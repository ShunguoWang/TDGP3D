function figs = tdgp3d_plotModel(P, varargin)
% TDGP3D_PLOTMODEL  Depth slices and vertical sections of a 3D TDGP posterior model.
%
%   P = tdgp3d_posteriorModels('.', 'burnin', 2000);
%   tdgp3d_plotModel(P)
%   tdgp3d_plotModel(P, 'depths', [500 1000 2000], 'x0', 0, 'y0', 0)
%   tdgp3d_plotModel('tdgp3d_out.posterior.mat')
%
%   P  the output of tdgp3d_posteriorModels, or its .mat file
%
% Upper panels: posterior mean of log10 resistivity. Lower panels: width of the 90%
% credible interval (95th minus 5th percentile, in decades); where it is large, the
% data do not constrain the model. MT sites are drawn as triangles and the domain
% for the GP node locations as a dashed line (outside it, the model is the
% background resistivity of MT3D_V6.dat).
%
% Two figures, saved as PNG in outDir:
%   <prefix>.depthSlices.png  maps at 'depths' (east y across, north x up)
%   <prefix>.sections.png     vertical sections along x at y = y0, and along y at x = x0
%
% Options (name/value):
%   'depths'  depths of the maps [m]           default: 4 depths between 10% and 70% of the domain
%   'x0'      x of the section along y [m]     default: median x of the sites
%   'y0'      y of the section along x [m]     default: median y of the sites
%   'zmax'    deepest depth of the sections [m]  default: bottom of the domain
%   'clim'    colour limits of the mean [log10 ohm-m]  default: range of the mean in the domain
%   'outDir', 'prefix'  default: as in tdgp3d_posteriorModels
%
% Output: figs = [figure of the depth slices, figure of the sections].
% Add the GP node locations with tdgp3d_plotNodes(runDir, 'last', figs).

if ischar(P)
    matFile = P;
    P = load(matFile);
    if isempty(P.outDir) || exist(P.outDir, 'dir') ~= 7
        P.outDir = fileparts(matFile);
    end
end

dom = P.domain;
if isempty(dom)
    dom = [min(P.x) max(P.x) min(P.y) max(P.y) 0 max(P.z)];     % unknown: the whole mesh
end
if isempty(P.sites.x)
    sx0 = mean(P.x);
    sy0 = mean(P.y);
else
    sx0 = median(P.sites.x);
    sy0 = median(P.sites.y);
end

p = inputParser;
p.FunctionName = 'tdgp3d_plotModel';
addParameter(p, 'depths', [],   @isnumeric);
addParameter(p, 'x0',     sx0,  @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'y0',     sy0,  @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'zmax',   dom(6), @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'clim',   [],   @(x) isnumeric(x) && (isempty(x) || numel(x) == 2));
addParameter(p, 'outDir', P.outDir, @ischar);
addParameter(p, 'prefix', P.prefix, @ischar);
parse(p, varargin{:});
o = p.Results;
if isempty(o.outDir)
    o.outDir = pwd;
end

depths = o.depths(:).';
if isempty(depths)
    depths = unique(round(dom(5) + (dom(6) - dom(5)) * [0.1 0.3 0.5 0.7], -1));
end
if any(depths < P.z(1) | depths > P.z(end))
    error('tdgp3d_plotModel:depth', 'depths must be between %g and %g m (the earth nodes of the mesh).', P.z(1), P.z(end));
end
zmax = min(o.zmax, P.z(end));

% --- mean and 90% interval width -----------------------------------------------
V  = P.log10rho_mean;
i5 = find(abs(P.percentiles - 5) < 1e-9, 1);
i95 = find(abs(P.percentiles - 95) < 1e-9, 1);
hasW = ~isempty(i5) && ~isempty(i95);
if hasW
    W = P.log10rho_pct(:,:,:,i95) - P.log10rho_pct(:,:,:,i5);
else
    warning('tdgp3d_plotModel:pct', 'The 5th and 95th percentiles are not in P: only the mean is plotted.');
end

if isempty(o.clim)
    [Zg, Xg, Yg] = ndgrid(P.z, P.x, P.y);
    in = Xg >= dom(1) & Xg <= dom(2) & Yg >= dom(3) & Yg <= dom(4) & Zg >= dom(5) & Zg <= dom(6);
    o.clim = [floor(10*min(V(in)))/10, ceil(10*max(V(in)))/10];
    if o.clim(2) <= o.clim(1)
        o.clim = o.clim + [-0.1 0.1];
    end
end
cmapM = flipud(jet(256));         % conductive = red, resistive = blue
cmapW = parula(256);
wlim  = [];
if hasW
    wlim = [0, max(0.1, ceil(10*max(W(:)))/10)];
end
nrow = 1 + hasW;
titleStr = sprintf('%s: posterior mean (top) and 90%% interval width p95-p05 (bottom), %d samples', ...
                   o.prefix, P.nSamples);
if ~hasW
    titleStr = sprintf('%s: posterior mean, %d samples', o.prefix, P.nSamples);
end

% --- depth slices -------------------------------------------------------------------
nd  = numel(depths);
fig1 = figure('Color', 'w', 'Position', [50 50 330*nd + 120, 300*nrow + 70], 'Visible', 'off');   % shown once saved
tl  = tiledlayout(fig1, nrow, nd, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:nrow
    for j = 1:nd
        ax = nexttile(tl);
        if r == 1
            A = local_zslice(V, P.z, depths(j));
        else
            A = local_zslice(W, P.z, depths(j));
        end
        pcolor(ax, P.y/1e3, P.x/1e3, A);                 % rows of A: x (north), columns: y (east)
        shading(ax, 'interp');
        hold(ax, 'on');
        plot(ax, dom([3 4 4 3 3])/1e3, dom([1 1 2 2 1])/1e3, 'k--', 'LineWidth', 0.8);
        plot(ax, P.sites.y/1e3, P.sites.x/1e3, 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
        axis(ax, 'equal', 'tight');
        title(ax, sprintf('depth %g m', depths(j)));
        set(ax, 'Tag', 'tdgp3dMap', 'UserData', depths(j));        % for tdgp3d_plotNodes
        local_style(ax, r == 2, cmapM, cmapW, o.clim, wlim);
        xlabel(ax, 'east y (km)');
        if j == 1
            ylabel(ax, 'north x (km)');
        end
        if j == nd
            local_colorbar(ax, r);
        end
    end
end
title(tl, titleStr, 'Interpreter', 'none');
tdgp3d_savePNG(fig1, fullfile(o.outDir, [o.prefix '.depthSlices.png']));
set(fig1, 'UserData', fullfile(o.outDir, [o.prefix '.depthSlices.png']));   % for tdgp3d_plotNodes

% --- vertical sections ----------------------------------------------------------------
span = max(P.x(end) - P.x(1), P.y(end) - P.y(1));
rowH = round(430 * min(max(1.04*zmax / span, 0.25), 1.5)) + 90;    % panel height from the section's aspect ratio
fig2 = figure('Color', 'w', 'Position', [80 80 1100, rowH*nrow + 70], 'Visible', 'off');   % shown once saved
tl  = tiledlayout(fig2, nrow, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:nrow
    if r == 1
        F = V;
    else
        F = W;
    end
    % along x at y = y0:
    ax = nexttile(tl);
    A = local_yslice(F, P.y, o.y0);                      % [nze x nx]
    pcolor(ax, P.x/1e3, P.z/1e3, A);
    shading(ax, 'interp');
    hold(ax, 'on');
    plot(ax, dom([1 2 2 1 1])/1e3, dom([5 5 6 6 5])/1e3, 'k--', 'LineWidth', 0.8);
    on = abs(P.sites.y - o.y0) <= 0.5*local_step(P.y);
    plot(ax, P.sites.x(on)/1e3, zeros(nnz(on), 1), 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal');
    xlim(ax, [P.x(1) P.x(end)]/1e3);
    ylim(ax, [-0.04*zmax zmax]/1e3);                      % a little room above the surface for the sites
    title(ax, sprintf('along x at y = %g km', o.y0/1e3));
    set(ax, 'Tag', 'tdgp3dSectionX', 'UserData', o.y0);            % for tdgp3d_plotNodes
    xlabel(ax, 'north x (km)');
    ylabel(ax, 'depth (km)');
    local_style(ax, r == 2, cmapM, cmapW, o.clim, wlim);
    % along y at x = x0:
    ax = nexttile(tl);
    A = local_xslice(F, P.x, o.x0);                      % [nze x ny]
    pcolor(ax, P.y/1e3, P.z/1e3, A);
    shading(ax, 'interp');
    hold(ax, 'on');
    plot(ax, dom([3 4 4 3 3])/1e3, dom([5 5 6 6 5])/1e3, 'k--', 'LineWidth', 0.8);
    on = abs(P.sites.x - o.x0) <= 0.5*local_step(P.x);
    plot(ax, P.sites.y(on)/1e3, zeros(nnz(on), 1), 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal');
    xlim(ax, [P.y(1) P.y(end)]/1e3);
    ylim(ax, [-0.04*zmax zmax]/1e3);                      % a little room above the surface for the sites
    title(ax, sprintf('along y at x = %g km', o.x0/1e3));
    set(ax, 'Tag', 'tdgp3dSectionY', 'UserData', o.x0);            % for tdgp3d_plotNodes
    xlabel(ax, 'east y (km)');
    local_style(ax, r == 2, cmapM, cmapW, o.clim, wlim);
    local_colorbar(ax, r);
end
title(tl, titleStr, 'Interpreter', 'none');
tdgp3d_savePNG(fig2, fullfile(o.outDir, [o.prefix '.sections.png']));
set(fig2, 'UserData', fullfile(o.outDir, [o.prefix '.sections.png']));      % for tdgp3d_plotNodes

fprintf('tdgp3d_plotModel: depth slices at %s m, sections at y = %g m and x = %g m, mean colour range %g to %g\n', ...
        mat2str(depths), o.y0, o.x0, o.clim(1), o.clim(2));
fprintf('  saved %s.depthSlices.png and %s.sections.png in %s\n', o.prefix, o.prefix, o.outDir);
figs = [fig1 fig2];

end


function local_style(ax, isWidth, cmapM, cmapW, clim, wlim)
% Colour map and limits: the mean (upper panels) or the interval width (lower panels).
if isWidth
    colormap(ax, cmapW);
    caxis(ax, wlim);      %#ok<CAXIS> -- clim() needs R2022a+; caxis works on older MATLAB too
else
    colormap(ax, cmapM);
    caxis(ax, clim);      %#ok<CAXIS>
end
set(ax, 'Layer', 'top', 'Box', 'on', 'TickDir', 'out');
end


function local_colorbar(ax, r)
cb = colorbar(ax);
if r == 1
    ylabel(cb, 'log_{10} resistivity (\Omega m)');
else
    ylabel(cb, 'p95 - p05 (decades)');
end
end


function A = local_zslice(F, z, d)
% [nx x ny] slice of F [nze x nx x ny] at depth d, linear in z between node levels.
i = find(z <= d, 1, 'last');
i = min(i, numel(z) - 1);
w = (d - z(i)) / (z(i+1) - z(i));
A = (1 - w) * reshape(F(i,:,:), size(F,2), size(F,3)) + w * reshape(F(i+1,:,:), size(F,2), size(F,3));
end


function A = local_yslice(F, y, y0)
% [nze x nx] section of F [nze x nx x ny] at y = y0, linear in y between node columns.
y0 = min(max(y0, y(1)), y(end));
j = min(find(y <= y0, 1, 'last'), numel(y) - 1);
w = (y0 - y(j)) / (y(j+1) - y(j));
A = (1 - w) * F(:,:,j) + w * F(:,:,j+1);
end


function A = local_xslice(F, x, x0)
% [nze x ny] section of F [nze x nx x ny] at x = x0, linear in x between node columns.
x0 = min(max(x0, x(1)), x(end));
i = min(find(x <= x0, 1, 'last'), numel(x) - 1);
w = (x0 - x(i)) / (x(i+1) - x(i));
A = (1 - w) * reshape(F(:,i,:), size(F,1), size(F,3)) + w * reshape(F(:,i+1,:), size(F,1), size(F,3));
end


function s = local_step(v)
if numel(v) > 1
    s = v(2) - v(1);
else
    s = Inf;
end
end
