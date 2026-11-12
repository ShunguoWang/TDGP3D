function P = tdgp3d_posteriorModels(runDir, varargin)
% TDGP3D_POSTERIORMODELS  Posterior mean and percentile models from a 3D TDGP ensemble.
%
%   P = tdgp3d_posteriorModels(runDir, 'burnin', 2000)
%   tdgp3d_plotModel(P)
%
%   runDir  folder with tdgp.settings, MT3D_V6.dat, the ModEM data file and
%           <outputBase>_ensemble.txt and _domain.txt (a copy of the cluster run folder)
%
% Every saved sample (a set of GP nodes) is evaluated at the SAME points, the earth
% nodes of the forward mesh (where the sampler evaluated it), so the samples can be
% compared node by node although the number of GP nodes k changes from sample to
% sample. At every earth node, statistics are taken across ALL the samples left after
% burn-in (and 'thin', if given): the mean of log10(rho) (i.e. the geometric mean of
% rho) and percentiles. They are computed for a block of mesh columns at a time, so the
% memory stays at about 1 GB however many samples there are. The result is saved in
% <outDir>/<prefix>.posterior.mat.
%
% Options (name/value):
%   'burnin'       discard samples with iteration <= burnin              [0]
%   'thin'         keep every thin-th sample after burn-in               [1]
%   'maxSamples'   use only the first maxSamples of them                 [Inf = all]
%   'percentiles'  percentiles, 0..100                                   [5 25 50 75 95]
%   'settingsFile' default <runDir>/tdgp.settings
%   'ensembleFile' default <runDir>/<outputBase>_ensemble.txt
%   'outDir'       default runDir
%   'prefix'       default outputBase from tdgp.settings
%
% Choose 'burnin' from <outputBase>_summary.txt (logL and k versus iteration).
%
% Output struct P:
%   P.x [nx x 1], P.y [ny x 1], P.z [nze x 1]  earth-node coordinates (m): x north, y east, z depth
%   P.log10rho_mean [nze x nx x ny]            posterior mean of log10(ohm-m)
%   P.percentiles [nP x 1], P.log10rho_pct [nze x nx x ny x nP]
%   P.nSamples, P.imcmc, P.k                   samples used, their iterations and node counts
%   P.sites  (fields code, x, y; empty if the data file was not found)
%   P.domain [xmin xmax ymin ymax zmin zmax]   domain for the GP node locations ([] if unknown)
%   P.logrho_bg                                GP prior mean = log10 background of MT3D_V6.dat
%   P.m0 [nSamples x 1]                        each sample's GP prior mean: its sampled background
%                                              (sample_m0 runs), otherwise P.logrho_bg
%   P.settings, P.outDir, P.prefix

p = inputParser;
p.FunctionName = 'tdgp3d_posteriorModels';
addParameter(p, 'burnin',       0,      @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'thin',         1,      @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'maxSamples',   Inf,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'percentiles',  [5 25 50 75 95], @(x) isnumeric(x) && ~isempty(x) && all(x >= 0 & x <= 100));
addParameter(p, 'settingsFile', '',     @ischar);
addParameter(p, 'ensembleFile', '',     @ischar);
addParameter(p, 'outDir',       '',     @ischar);
addParameter(p, 'prefix',       '',     @ischar);
parse(p, varargin{:});
o = p.Results;

[S, F] = tdgp3d_resolveRun(runDir, o);
if exist(F.outDir, 'dir') ~= 7
    mkdir(F.outDir);
end

M = tdgp3d_readMesh(F.mesh);
if isfield(S, 'logrho_bg') && ~isempty(F.domain) && abs(S.logrho_bg - M.logrho_bg) > 1e-6
    warning('tdgp3d_posteriorModels:background', ...
            ['The run used a background of log10(rho) = %g (%s) but %s has %g: ' ...
             'is this the MT3D_V6.dat of the run?'], S.logrho_bg, F.domain, F.mesh, M.logrho_bg);
end
lambda = [S.lambda_x S.lambda_y S.lambda_z];

fprintf('tdgp3d_posteriorModels\n');
fprintf('  mesh      %s  (%d x %d x %d earth nodes)\n', F.mesh, numel(M.x), numel(M.y), numel(M.z) - M.jf + 1);
fprintf('  ensemble  %s\n', F.ensemble);
fprintf('  GP        lambda_x = %g, lambda_y = %g, lambda_z = %g m, sigma_m = %g, prior mean log10(rho) = %g\n', ...
        lambda, S.sigma_m, M.logrho_bg);

% --- the earth nodes, in the order the sampler used ([nze x nx x ny], z fastest) ---
% They form a tensor grid, so each axis is kept separately, in length scales; see
% local_gridMean.
ze  = M.z(M.jf:end);
nze = numel(ze);  nx = numel(M.x);  ny = numel(M.y);
xg  = M.x(:) / lambda(1);
yg  = M.y(:) / lambda(2);
zg  = ze(:) / lambda(3);
inc = S.lambda_increase;
if inc > 0
    fprintf('  lambda_increase = %g: length scales %g, %g, %g m at the surface, %g, %g, %g m at %g m\n', ...
            inc, lambda, lambda + inc*max(ze(end), 0), ze(end));
end

E = tdgp3d_readEnsemble(F.ensemble, 'burnin', o.burnin, 'thin', o.thin, 'maxSamples', o.maxSamples);
n = numel(E.imcmc);
if n == 0
    error('tdgp3d_posteriorModels:noSamples', ...
          'No samples left after burn-in = %g. Check ''burnin'' against the iterations in %s', o.burnin, F.ensemble);
end
fprintf('  %d samples read (iterations %d to %d, burn-in %g, thin %d)\n', n, E.imcmc(1), E.imcmc(end), o.burnin, o.thin);
m0s = E.m0;                                          % each sample's background (sample_m0 runs) ...
m0s(isnan(m0s)) = M.logrho_bg;                       % ... or that of MT3D_V6.dat
if any(~isnan(E.m0))
    fprintf('  sampled background m0: log10(rho) %.4f to %.4f (median %.4f)\n', min(m0s), max(m0s), median(m0s));
end

% --- GP weights of every sample (one Cholesky factorisation each) ----------------
nodeS = cell(n, 1);                                  % node positions
alpha = cell(n, 1);
good  = false(n, 1);
for i = 1:n
    nodes = E.nodes{i};
    [~, good(i), alpha{i}] = tdgp3d_gpMean(nodes, [], [], [], lambda, S.sigma_m, m0s(i), inc);
    if inc > 0
        nodeS{i} = nodes(:,1:3);                     % in m (local_gridMeanInc)
    else
        nodeS{i} = nodes(:,1:3) ./ lambda;           % in length scales (local_gridMean)
    end
end
if ~all(good)
    warning('tdgp3d_posteriorModels:gp', ...
            '%d of %d samples had a non-positive-definite kernel matrix and were dropped.', nnz(~good), n);
end
nodeS = nodeS(good);
alpha = alpha(good);
m0g   = m0s(good);
ng    = nnz(good);
if ng == 0
    error('tdgp3d_posteriorModels:noSamples', 'No usable samples.');
end

% --- statistics across the samples, a block of mesh columns at a time -------------
% All samples at all earth nodes would take ng*nze*nx*ny*4 bytes (10,000 samples on a
% 42 x 42 x 45 mesh: 3.2 GB, and as much again to sort them); a block takes ~0.5 GB.
pct  = o.percentiles(:);
ncol = nx * ny;
nper = max(1, min(ncol, floor(0.5e9 / (4 * ng * nze))));     % mesh columns per block
nblk = ceil(ncol / nper);
avg  = zeros(nze * ncol, 1);
Q    = zeros(nze * ncol, numel(pct));
for b = 1:nblk
    cols = (b-1)*nper + 1 : min(ncol, b*nper);
    ix   = mod(cols(:) - 1, nx) + 1;
    iy   = floor((cols(:) - 1) / nx) + 1;
    Xb   = zeros(nze * numel(cols), ng, 'single');            % one column per sample
    for i = 1:ng
        if inc > 0
            V = local_gridMeanInc(nodeS{i}, alpha{i}, M.x, M.y, ze, lambda, inc, ix, iy, m0g(i));
        else
            V = local_gridMean(nodeS{i}, alpha{i}, xg, yg, zg, ix, iy, m0g(i));
        end
        Xb(:, i) = V(:);
    end
    lin       = (cols(1) - 1) * nze + 1 : cols(end) * nze;    % z fastest, then x, then y
    avg(lin)  = mean(Xb, 2, 'double');
    Q(lin, :) = local_percentile(Xb, pct);
    if nblk > 1
        fprintf('  mesh columns %d-%d of %d done\n', cols(1), cols(end), ncol);
    end
end

sz = [nze nx ny];
P.x             = M.x;
P.y             = M.y;
P.z             = ze;
P.log10rho_mean = reshape(avg, sz);
P.percentiles   = pct;
P.log10rho_pct  = reshape(Q, [sz numel(pct)]);
P.nSamples      = ng;
P.imcmc         = E.imcmc(good);
P.k             = E.k(good);
P.sites         = struct('code', {{}}, 'x', zeros(0,1), 'y', zeros(0,1));
if ~isempty(F.data)
    P.sites = tdgp3d_readSites(F.data);
end
P.domain = [];
if ~isempty(F.domain)
    P.domain = [S.xmin S.xmax S.ymin S.ymax S.zmin S.zmax];
end
P.logrho_bg = M.logrho_bg;
P.m0        = m0g;
P.settings  = S;
P.outDir    = F.outDir;
P.prefix    = F.prefix;

matFile = fullfile(F.outDir, [F.prefix '.posterior.mat']);
save(matFile, '-struct', 'P');
fprintf('  used %d samples, k = %d to %d nodes (median %g)\n', P.nSamples, min(P.k), max(P.k), median(P.k));
fprintf('  saved %s\n', matFile);

end


function V = local_gridMean(nodeS, alpha, xg, yg, zg, ix, iy, mean0)
% GP mean of one sample at the earth nodes of the mesh columns (xg(ix), yg(iy)):
% V [nze x ncolumns]. nodeS [k x 3] node positions and xg, yg, zg the mesh axes, all in
% length scales; alpha [k x 1] the GP weights from tdgp3d_gpMean.
% The kernel exp(-0.5*(dx^2 + dy^2 + dz^2)) is a product of 1D Gaussians, so the value
% tdgp3d_gpMean gives, mean0 + sum_j alpha_j K(point, node j), is here
%   mean0 + sum_j alpha_j Kz(z, j) Kx(x, j) Ky(y, j)   (one matrix product).
Kx = exp(-0.5 * (xg - nodeS(:,1).').^2);           % [nx x k]
Ky = exp(-0.5 * (yg - nodeS(:,2).').^2);           % [ny x k]
Kz = exp(-0.5 * (zg - nodeS(:,3).').^2);           % [nze x k]
V  = mean0 + (Kz .* alpha.') * (Kx(ix, :) .* Ky(iy, :)).';
end


function V = local_gridMeanInc(nodes, alpha, xm, ym, ze, lambda, inc, ix, iy, mean0)
% local_gridMean for runs with lambda_increase > 0, every length scale lambda + inc*depth
% (tdgp3d_gpMean). nodes [k x 3] and the mesh axes xm, ym, ze in m. At one depth level
% the kernel is still a product of x, y and z factors, so each level is one matrix
% product (pagemtimes over the levels).
k   = size(nodes, 1);
nze = numel(ze);
[uy, ~, jy] = unique(iy);                                          % the y lines of this block
ln  = lambda(:).' + inc * max(nodes(:,3), 0);                      % [k x 3] at the nodes
lg  = lambda(:).' + inc * max(ze(:), 0);                           % [nze x 3] at the depth levels
lgx = reshape(lg(:,1), 1, 1, nze);  lgy = reshape(lg(:,2), 1, 1, nze);
sx  = ln(:,1).'.^2 + lgx.^2;                                       % [1 x k x nze]
sy  = ln(:,2).'.^2 + lgy.^2;
sz  = ln(:,3).'.^2 + lg(:,3).^2;                                   % [nze x k]
Kx  = sqrt(2 * ln(:,1).' .* lgx ./ sx) .* exp(-(xm(:) - nodes(:,1).').^2 ./ sx);        % [nx x k x nze]
Ky  = sqrt(2 * ln(:,2).' .* lgy ./ sy) .* exp(-(ym(uy(:)) - nodes(:,2).').^2 ./ sy);    % [nyb x k x nze]
Kz  = sqrt(2 * ln(:,3).' .* lg(:,3) ./ sz) .* exp(-(ze(:) - nodes(:,3).').^2 ./ sz);    % [nze x k]
A   = reshape((Kz .* alpha.').', 1, k, nze);
C   = permute(pagemtimes(Kx, permute(Ky .* A, [2 1 3])), [3 1 2]);                     % [nze x nx x nyb]
V   = mean0 + C(:, sub2ind([numel(xm) numel(uy)], ix(:), jy(:)));
end


function Q = local_percentile(X, p)
% Percentiles of the rows of X (columns = samples): linear interpolation between
% order statistics, position 1 + (n-1)*p/100 (as in the 2D tdgp_posteriorModels).
n  = size(X, 2);
Xs = sort(X, 2);
Q  = zeros(size(X, 1), numel(p));
for i = 1:numel(p)
    pos = 1 + (n - 1) * p(i) / 100;
    lo  = floor(pos);
    hi  = ceil(pos);
    w   = pos - lo;
    Q(:, i) = (1 - w) * double(Xs(:, lo)) + w * double(Xs(:, hi));
end
end
