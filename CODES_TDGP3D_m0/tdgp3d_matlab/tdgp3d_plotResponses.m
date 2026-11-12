function out = tdgp3d_plotResponses(runDir, varargin)
% TDGP3D_PLOTRESPONSES  Observed data and the predicted responses of saved samples, along profiles.
%
%   out = tdgp3d_plotResponses('.', 'burnin', 300)
%   out = tdgp3d_plotResponses('.', 'burnin', 300, 'x0', 100, 'y0', -100)
%
%   runDir  folder with tdgp.settings, the ModEM data file and <outputBase>_responses.txt
%
% The sampler saves the predicted data of every 10th iteration of the beta = 1 chain in
% <outputBase>_responses.txt. For the saved iterations after 'burnin' this plots, along two
% profiles of sites -- along x at y = y0, and along y at x = x0 -- the apparent resistivity and
% phase of Zxy and Zyx and the tipper Tzx and Tzy (real and imaginary parts), one column per
% frequency: the observed data with their error bars (black), the responses of the samples (grey)
% and the sample with the smallest misfit (red). Phases and the tipper are shown for time dependence
% exp(+i omega t), i.e. Zxy near +45 degrees; 180 degrees is added to the phase of Zyx.
% Saved as <prefix>.responses_x.png and <prefix>.responses_y.png in outDir.
%
% Options (name/value):
%   'burnin'      leave out iterations <= burnin                            [0]
%   'thin'        use every thin-th saved response                          [1]
%   'maxSamples'  use only the first maxSamples of them                     [Inf = all]
%   'x0'          x of the profile along y [m]: the column of sites nearest to it
%   'y0'          y of the profile along x [m]: the row of sites nearest to it   [middle of the sites]
%   'outDir', 'prefix'  default: runDir and outputBase from tdgp.settings
%
% Output struct: out.imcmc (iterations used), out.rms (RMS misfit of each), out.best (iteration
% with the smallest misfit), out.rmsTable (RMS of that sample per component and frequency, a table),
% out.files.

p = inputParser;
p.FunctionName = 'tdgp3d_plotResponses';
addParameter(p, 'burnin',     0,   @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'thin',       1,   @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'maxSamples', Inf, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'x0',         [],  @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addParameter(p, 'y0',         [],  @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addParameter(p, 'outDir',     '',  @ischar);
addParameter(p, 'prefix',     '',  @ischar);
parse(p, varargin{:});
o = p.Results;

[S, F] = tdgp3d_resolveRun(runDir, struct('settingsFile', '', 'ensembleFile', '', 'outDir', o.outDir, 'prefix', o.prefix));
if isempty(F.data)
    error('tdgp3d_plotResponses:noData', 'Data file %s not found in %s', S.datafile, runDir);
end
if exist(F.outDir, 'dir') ~= 7
    mkdir(F.outDir);
end
D  = tdgp3d_readModEM(F.data);
nr = numel(D.period);
nd = 2 * nr;
dobs = reshape([D.re D.im].', [], 1);
dunc = reshape([D.err D.err].', [], 1);

% --- the saved responses ---------------------------------------------------------------
respFile = fullfile(runDir, [S.outputbase '_responses.txt']);
fid = fopen(respFile, 'r');
if fid < 0
    error('tdgp3d_plotResponses:noResponses', 'Cannot open %s', respFile);
end
it = zeros(256, 1);  R = zeros(256, nd);  m = 0;  nElig = 0;  nBad = 0;  hitMax = false;
while true
    L = fgetl(fid);
    if ~ischar(L), break; end
    v = sscanf(L, '%f');
    if numel(v) ~= nd + 1
        if ~isempty(v), nBad = nBad + 1; end          % a line cut short, or another data file
        continue;
    end
    if v(1) <= o.burnin, continue; end
    nElig = nElig + 1;
    if mod(nElig - 1, o.thin) ~= 0, continue; end
    if m >= o.maxSamples, hitMax = true; break; end
    m = m + 1;
    if m > numel(it)                                   % grow by doubling
        it(2*m, 1)  = 0;
        R(2*m, nd)  = 0;
    end
    it(m)    = v(1);
    R(m, :)  = v(2:end).';
end
fclose(fid);
it = it(1:m);
R  = R(1:m, :);
if isempty(it)
    error('tdgp3d_plotResponses:noSamples', ...
          'No responses after burn-in = %g in %s (%d lines do not match the %d data of %s).', ...
          o.burnin, respFile, nBad, nd, F.data);
end
if nBad > 0
    warning('tdgp3d_plotResponses:badLines', '%d line(s) of %s skipped (cut short, or not %d values).', nBad, respFile, nd);
end
if hitMax
    warning('tdgp3d_plotResponses:maxSamples', 'Only the first %d samples are drawn; raise ''thin'' to cover the run.', o.maxSamples);
end
chi2 = sum(((R - dobs.') ./ dunc.').^2, 2);
rms  = sqrt(chi2 / nd);
[~, ib] = min(chi2);
fprintf('tdgp3d_plotResponses: %d samples (iterations %d to %d); smallest misfit at iteration %d, RMS %.3f\n', ...
        numel(it), it(1), it(end), it(ib), rms(ib));

% --- RMS of the best sample per component and frequency ----------------------------------------
per   = unique(D.period);
comps = {'ZXX', 'ZXY', 'ZYX', 'ZYY', 'TX', 'TY'};
T = nan(numel(comps), numel(per));
res = (R(ib, :).' - dobs) ./ dunc;
for c = 1:numel(comps)
    for k = 1:numel(per)
        rr = find(strcmp(D.comp, comps{c}) & abs(D.period - per(k)) <= 1e-6 * per(k));
        if ~isempty(rr)
            T(c, k) = sqrt(mean(res([2*rr-1; 2*rr]).^2));
        end
    end
end
rmsTable = array2table(T, 'RowNames', comps, 'VariableNames', strrep(compose('f_%gHz', 1 ./ per.'), '.', 'p'));
fprintf('  RMS of iteration %d per component (rows) and frequency (columns):\n', it(ib));
disp(rmsTable);

% --- profiles ------------------------------------------------------------------------------
[codes, iu] = unique(D.code, 'stable');
sx = D.x(iu);
sy = D.y(iu);
if isempty(o.x0), o.x0 = median(sx); end
if isempty(o.y0), o.y0 = median(sy); end
ux = unique(sx);  [~, i] = min(abs(ux - o.x0));  x0row = ux(i);
uy = unique(sy);  [~, i] = min(abs(uy - o.y0));  y0row = uy(i);

C = struct('D', D, 'R', R, 'it', it, 'ib', ib, 'rms', rms, 'per', per, 'codes', {codes}, 'F', F);
files = {};
files{end+1} = local_profile(C, 'x', abs(sy - y0row) < 0.5, sx, y0row);
files{end+1} = local_profile(C, 'y', abs(sx - x0row) < 0.5, sy, x0row);

out.imcmc    = it;
out.rms      = rms;
out.best     = it(ib);
out.rmsTable = rmsTable;
out.files    = files;

end


function f = local_profile(C, dirn, onLine, spos, other)
% One figure: the sites of onLine along direction dirn ('x' or 'y'; spos = their coordinate),
% for the data C.D, the saved responses C.R (iterations C.it, best C.ib) and the periods C.per.
D = C.D;  R = C.R;  it = C.it;  ib = C.ib;  rms = C.rms;  per = C.per;  codes = C.codes;  F = C.F;
    idx   = find(onLine);
    [pos, o2] = sort(spos(idx));
    idx   = idx(o2);
    pos   = pos / 1e3;
    nf    = numel(per);
    mu0   = 4e-7 * pi;
    lab   = {'\rho_a xy (\Omega m)', '\phi xy (deg)', '\rho_a yx (\Omega m)', '\phi yx + 180 (deg)', 'T_{zx}', 'T_{zy}'};
    cmpOf = {'ZXY', 'ZXY', 'ZYX', 'ZYX', 'TX', 'TY'};
    fig = figure('Color', 'w', 'Position', [50 50 320*nf + 120, 980], 'Visible', 'off');   % shown once saved
    tl  = tiledlayout(fig, 6, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
    for q = 1:6
        for k = 1:nf
            ax = nexttile(tl, (q-1)*nf + k);
            hold(ax, 'on');
            obs = nan(numel(idx), 1); ebl = obs; ebu = obs; P = nan(numel(it), numel(idx));
            obsI = obs; PI = P;                                   % imaginary parts (tipper)
            for j = 1:numel(idx)
                rr = find(strcmp(D.code, codes{idx(j)}) & strcmp(D.comp, cmpOf{q}) & ...
                          abs(D.period - per(k)) <= 1e-6 * per(k), 1);
                if isempty(rr), continue; end
                zo = D.re(rr) + 1i * D.sgn(rr) * D.im(rr);                  % exp(+i omega t)
                zp = R(:, 2*rr-1) + 1i * D.sgn(rr) * R(:, 2*rr);
                if q <= 4
                    a  = abs(zo) * D.toOhm(rr);  da = D.err(rr) * D.toOhm(rr);
                    if mod(q, 2) == 1                                         % apparent resistivity
                        obs(j) = a^2 * per(k) / (2*pi*mu0);
                        ebu(j) = 2 * obs(j) * da / a;
                        ebl(j) = min(ebu(j), 0.9 * obs(j));                  % keep the bar above 0 on the log axis
                        P(:, j) = (abs(zp) * D.toOhm(rr)).^2 * per(k) / (2*pi*mu0);
                    else                                                      % phase
                        add = 180 * strcmp(cmpOf{q}, 'ZYX');
                        obs(j) = local_wrap(angle(zo) * 180/pi + add);
                        ebu(j) = min(da / a, 1) * 180/pi;  ebl(j) = ebu(j);
                        P(:, j) = local_wrap(angle(zp) * 180/pi + add);
                    end
                else                                                          % tipper: Re and Im
                    obs(j) = real(zo);  obsI(j) = imag(zo);  ebu(j) = D.err(rr);  ebl(j) = ebu(j);
                    P(:, j) = real(zp); PI(:, j) = imag(zp);
                end
            end
            grey = [0.72 0.72 0.72];
            plot(ax, pos, P.', '-', 'Color', grey, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            if q >= 5
                plot(ax, pos, PI.', '--', 'Color', grey, 'LineWidth', 0.5);
                plot(ax, pos, P(ib, :), 'r-', 'LineWidth', 1.2);
                plot(ax, pos, PI(ib, :), 'r--', 'LineWidth', 1.2);
                errorbar(ax, pos, obs, ebl, ebu, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'LineStyle', 'none');
                errorbar(ax, pos, obsI, ebl, ebu, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineStyle', 'none');
            else
                hb = plot(ax, pos, P(ib, :), 'r-', 'LineWidth', 1.2);
                ho = errorbar(ax, pos, obs, ebl, ebu, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'LineStyle', 'none');
                if mod(q, 2) == 1
                    set(ax, 'YScale', 'log');
                end
                if k == 1 && q == 1                              % one legend for all panels, above them
                    hs = plot(ax, NaN, NaN, '-', 'Color', grey);
                    hi = plot(ax, NaN, NaN, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
                    hd = plot(ax, NaN, NaN, 'r--', 'LineWidth', 1.2);
                    lg = legend(ax, [ho hs hb hi hd], {'observed', 'samples', sprintf('best (it. %d)', it(ib)), ...
                                'T imaginary part: observed', 'T imaginary part: best'}, ...
                                'Orientation', 'horizontal', 'FontSize', 8);
                    lg.Layout.Tile = 'north';
                end
            end
            set(ax, 'Box', 'on', 'TickDir', 'out');
            xlim(ax, [min(pos) max(pos)] + [-0.05 0.05] * max(0.2, max(pos) - min(pos)));
            if q == 1, title(ax, sprintf('%g Hz', 1/per(k))); end
            if k == 1, ylabel(ax, lab{q}); end
            if q == 6
                if dirn == 'x', xlabel(ax, 'north x (km)'); else, xlabel(ax, 'east y (km)'); end
            end
        end
    end
    if dirn == 'x'
        where = sprintf('along x at y = %g km', other/1e3);
    else
        where = sprintf('along y at x = %g km', other/1e3);
    end
    title(tl, sprintf('%s: data and responses %s (%d sites); %d samples, iterations %d-%d; best RMS %.2f', ...
          F.prefix, where, numel(idx), numel(it), it(1), it(end), rms(ib)), 'Interpreter', 'none');
    f = fullfile(F.outDir, sprintf('%s.responses_%s.png', F.prefix, dirn));
    tdgp3d_savePNG(fig, f);
    fprintf('  saved %s\n', f);
end


function a = local_wrap(a)
% Angles in degrees to (-180, 180].
a = mod(a + 180, 360) - 180;
a(a == -180) = 180;
end
