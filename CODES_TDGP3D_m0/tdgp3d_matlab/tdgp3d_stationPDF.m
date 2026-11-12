function out = tdgp3d_stationPDF(runDir, varargin)
% TDGP3D_STATIONPDF  Posterior PDF of log10 resistivity versus depth beneath each MT site.
%
%   out = tdgp3d_stationPDF('.', 'burnin', 3000)
%   out = tdgp3d_stationPDF('.', 'burnin', 3000, 'sites', {'S0505', 'S0506', 'S0605', 'S0606'})
%
%   runDir  folder with tdgp.settings, MT3D_V6.dat, the ModEM data file and
%           <outputBase>_ensemble.txt (and <outputBase>_domain.txt)
%
% The 3D version of tdgp_stationPDF (2D code). Every sample left after burn-in (and
% 'thin', if given) is evaluated on a vertical profile beneath each site, from the surface down 'depth'
% metres. At each depth the values of the samples go into a histogram of log10
% resistivity (0.1 decade bins over the prior range), normalised to a PDF that
% integrates to 1 over log10 resistivity. One panel per site is plotted, with the
% median (solid) and the 5 and 95 percentiles (dashed), all panels sharing one PDF
% color scale (1/decade) and colorbar so they are directly comparable. When the sites
% form a full grid, the panels are placed as the sites are on the map (north at the
% top, east to the right). Saved as <prefix>.stationPDF.png; everything is saved in
% <prefix>.stationPDF.mat.
%
% Options (name/value):
%   'burnin'   discard samples with iteration <= burnin                  [0]
%   'thin'     keep every thin-th sample after burn-in                  [1]
%   'depth'    length of the profile below each site [m]  [bottom of the domain]
%   'sites'    site codes to plot, e.g. {'S0505', 'S0506'}        [all sites]
%   'outDir', 'prefix'   default: runDir and outputBase from tdgp.settings
%
% Output struct (ns sites, nz depths, nb bins):
%   out.code {ns x 1}, out.x, out.y [ns x 1]   site codes and positions (m; x north, y east)
%   out.depth  [nz x 1]          depth (m); below the domain the PDF is NaN
%   out.logrho [nb x 1]          bin centres, log10(ohm-m)
%   out.pdf    [nz x nb x ns]    PDF (1/decade)
%   out.median, out.p05, out.p95 [nz x ns]   log10(ohm-m)
%   out.nSamples                 samples used

p = inputParser;
p.FunctionName = 'tdgp3d_stationPDF';
addParameter(p, 'burnin', 0,  @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'thin',   1,  @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'depth',  [], @(x) isnumeric(x) && (isempty(x) || (isscalar(x) && x > 0)));
addParameter(p, 'sites',  {}, @(x) iscellstr(x) || ischar(x) || isstring(x));
addParameter(p, 'outDir', '', @ischar);
addParameter(p, 'prefix', '', @ischar);
parse(p, varargin{:});
o = p.Results;

[S, F] = tdgp3d_resolveRun(runDir, struct('settingsFile', '', 'ensembleFile', '', 'outDir', o.outDir, 'prefix', o.prefix));
if isempty(F.data)
    error('tdgp3d_stationPDF:noData', 'Data file %s not found in %s (the sites come from it)', S.datafile, runDir);
end
if exist(F.outDir, 'dir') ~= 7
    mkdir(F.outDir);
end

% --- GP prior mean and the bottom of the domain for the node locations ---------------
M     = tdgp3d_readMesh(F.mesh);
mean0 = M.logrho_bg;
zbot  = M.z(end);
if isfield(S, 'logrho_bg')
    mean0 = S.logrho_bg;                              % what the run used (<outputBase>_domain.txt)
end
if ~isempty(F.domain)
    zbot = S.zmax;
end
if isempty(o.depth)
    o.depth = zbot;
end
lambda = [S.lambda_x S.lambda_y S.lambda_z];

% --- sites, from the data file ---------------------------------------------------------
R = tdgp3d_readSites(F.data);
if isempty(o.sites)
    sel = (1:numel(R.code)).';
else
    want = cellstr(o.sites);
    [found, sel] = ismember(want(:), R.code);
    if ~all(found)
        error('tdgp3d_stationPDF:site', 'Site(s) %s not in %s', strjoin(want(~found).', ', '), F.data);
    end
end
xs = R.x(sel);  ys = R.y(sel);  code = R.code(sel);
ns = numel(sel);
ux = unique(xs);  uy = unique(ys);
isMap = numel(ux) > 1 && numel(uy) > 1 && numel(ux) * numel(uy) == ns && size(unique([xs ys], 'rows'), 1) == ns;
if isMap                                              % a full grid: panels as on the map
    [~, ord] = sortrows([-xs ys]);
    nrow = numel(ux);  ncol = numel(uy);
else
    ord  = (1:ns).';
    ncol = min(5, ns);  nrow = ceil(ns / ncol);
end
xs = xs(ord);  ys = ys(ord);  code = code(ord);

% --- profiles beneath the sites (on the surface, z = 0) --------------------------------
zd    = linspace(0, o.depth, 201).';
nz    = numel(zd);
Xpts  = repmat(xs.', nz, 1);
Ypts  = repmat(ys.', nz, 1);
Zpts  = repmat(zd, 1, ns);
valid = Zpts <= zbot + 1e-6;                          % inside the domain for the node locations

% --- histogram of log10 resistivity at every profile point ------------------------------
db    = 0.1;
nb    = round((S.logrho_max - S.logrho_min) / db);
edges = S.logrho_min + db * (0:nb);
ctr   = edges(1:end-1) + db/2;

E = tdgp3d_readEnsemble(F.ensemble, 'burnin', o.burnin, 'thin', o.thin);     % all samples after burn-in
n = numel(E.imcmc);
if n == 0
    error('tdgp3d_stationPDF:noSamples', 'No samples left after burn-in = %g in %s', o.burnin, F.ensemble);
end
fprintf('tdgp3d_stationPDF: %d samples (iterations %d to %d), %d sites, %d depths to %g m\n', ...
        n, E.imcmc(1), E.imcmc(end), ns, nz, o.depth);

m0s = E.m0;                                           % each sample's background (sample_m0 runs) ...
m0s(isnan(m0s)) = mean0;                              % ... or the run's fixed one
counts = zeros(nz, nb, ns);
[Zi, Ri] = ndgrid(1:nz, 1:ns);
ngood = 0;
for i = 1:n
    [mu, ok] = tdgp3d_gpMean(E.nodes{i}, Xpts(:), Ypts(:), Zpts(:), lambda, S.sigma_m, m0s(i), S.lambda_increase);
    if ~ok
        continue;
    end
    ngood = ngood + 1;
    ib  = min(max(floor((reshape(mu, nz, ns) - S.logrho_min) / db) + 1, 1), nb);
    lin = sub2ind([nz nb ns], Zi(:), ib(:), Ri(:));
    counts(lin) = counts(lin) + 1;                    % one count per depth and site, so no repeats
end
if ngood == 0
    error('tdgp3d_stationPDF:noSamples', 'No usable samples.');
end

pdf = counts / (ngood * db);
cdf = cumsum(counts, 2) / ngood;
for k = 1:ns
    pdf(~valid(:,k), :, k) = NaN;
end

out.code     = code;
out.x        = xs;
out.y        = ys;
out.depth    = zd;
out.logrho   = ctr(:);
out.pdf      = pdf;
out.median   = local_quantile(cdf, edges, 0.50, valid);
out.p05      = local_quantile(cdf, edges, 0.05, valid);
out.p95      = local_quantile(cdf, edges, 0.95, valid);
out.nSamples = ngood;

% --- plot and save ---------------------------------------------------------------------
pdfMax = max(pdf(:), [], 'omitnan');    % common color scale so every panel is directly comparable
fig = figure('Color', 'w', 'Position', [50 50 max(700, 150*ncol + 170) 140*nrow + 150], 'Visible', 'off');   % shown once saved
tl  = tiledlayout(fig, nrow, ncol, 'TileSpacing', 'tight', 'Padding', 'compact');
for k = 1:ns
    ax = nexttile(tl);
    imagesc(ax, ctr, zd/1e3, pdf(:,:,k)); hold(ax, 'on');
    colormap(ax, flipud(gray));
    caxis(ax, [0 pdfMax]);    %#ok<CAXIS> -- clim() needs R2022a+; caxis works on older MATLAB too
    hmed = plot(ax, out.median(:,k), zd/1e3, 'r-',  'LineWidth', 1.0);
    hci  = plot(ax, out.p05(:,k),    zd/1e3, 'r--', 'LineWidth', 0.6);
    plot(ax, out.p95(:,k), zd/1e3, 'r--', 'LineWidth', 0.6);
    title(ax, sprintf('%s: x %g, y %g km', code{k}, xs(k)/1e3, ys(k)/1e3), ...
          'FontSize', 7, 'FontWeight', 'normal', 'Interpreter', 'none');
    set(ax, 'YDir', 'reverse', 'XLim', [S.logrho_min S.logrho_max], 'FontSize', 7, ...
        'XTick', ceil(S.logrho_min + 0.5):floor(S.logrho_max - 0.5));   % no labels at the edges, where panels meet
    if k <= ns - ncol                                 % a panel below: no x tick labels
        set(ax, 'XTickLabel', []);
    end
    if mod(k-1, ncol) ~= 0
        set(ax, 'YTickLabel', []);
    end
    if k == 1
        lg = legend(ax, [hmed hci], {'median', '5th/95th percentile'}, 'Orientation', 'horizontal', 'FontSize', 8);
    end
end
lg.Layout.Tile = 'north';
cb = colorbar(ax);              % one shared colorbar for the whole layout (same scale on every panel)
cb.Layout.Tile = 'east';
ylabel(cb, 'PDF (1/decade)');
xlabel(tl, 'log_{10} resistivity (\Omega m)');
ylabel(tl, 'depth (km)');
where = '';
if isMap
    where = ' (panels placed as the sites are on the map: north at the top, east to the right)';
end
title(tl, sprintf('%s: posterior PDF beneath each site%s; %d samples, iterations %d-%d', ...
      F.prefix, where, ngood, E.imcmc(1), E.imcmc(end)), 'Interpreter', 'none', 'FontSize', 10);
tdgp3d_savePNG(fig, fullfile(F.outDir, [F.prefix '.stationPDF.png']));
save(fullfile(F.outDir, [F.prefix '.stationPDF.mat']), '-struct', 'out');
fprintf('  saved %s.stationPDF.png and .mat in %s\n', F.prefix, F.outDir);

end


function q = local_quantile(cdf, edges, pp, valid)
% pp-quantile of each histogram (cumulative fractions cdf along dimension 2),
% linear inside the bin; NaN where the profile point is outside the domain.
[nz, nb, ns] = size(cdf);
db = edges(2) - edges(1);
q  = nan(nz, ns);
for k = 1:ns
    for j = 1:nz
        if ~valid(j,k)
            continue;
        end
        c  = squeeze(cdf(j,:,k));
        m  = find(c >= pp, 1);
        if isempty(m)
            m = nb;
        end
        c0 = 0;
        if m > 1
            c0 = c(m-1);
        end
        q(j,k) = edges(m) + db * (pp - c0) / max(c(m) - c0, eps);
    end
end
end
