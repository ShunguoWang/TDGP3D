function h = tdgp3d_plotNodes(runDir, sel, figs)
% TDGP3D_PLOTNODES  Overlay GP node locations on the figures of tdgp3d_plotModel.
%
%   figs = tdgp3d_plotModel(P);
%   tdgp3d_plotNodes('.', 'last', figs)             % nodes of the last saved sample
%   tdgp3d_plotNodes('.', 'best', figs)             % nodes of the sample with the smallest misfit
%   tdgp3d_plotNodes('.', 700, figs)                % nodes of one iteration
%   tdgp3d_plotNodes('.', 400:735, figs)            % nodes of several samples
%
%   runDir  folder with tdgp.settings, MT3D_V6.dat and <outputBase>_ensemble.txt
%   sel     'last', 'best' (the iteration with the smallest chi-squared in
%           <outputBase>_summary.txt), or one or more iteration numbers (first
%           column of the ensemble file)
%   figs    the figures returned by tdgp3d_plotModel (default: all open figures)
%
% A node is a point in 3D, so each panel shows the nodes near what it shows: a map
% at depth d the nodes within one vertical length scale (lambda_z) of d, the section
% along x at y = y0 the nodes within lambda_y of y0, and the section along y at
% x = x0 those within lambda_x of x0. Farther nodes hardly affect the panel (GP
% kernel below exp(-0.5) = 0.61). With lambda_increase > 0 each window is the node's
% own length scale at its depth, lambda + lambda_increase*depth.
%
% One sample is drawn as white circles with black edges, several samples as small
% black dots (which shows where the nodes concentrate), as in the 2D tdgp_plotNodes.
% The figures are saved again with the nodes, as <name>.nodes.png next to the PNGs of
% tdgp3d_plotModel. The overlay carries the tag 'tdgpNodes': calling the function
% again replaces it, and  delete(findobj(figs, 'tag', 'tdgpNodes'))  removes it.
%
% Output h: handles of the plotted node markers.

if nargin < 3 || isempty(figs)
    figs = findobj(groot, 'Type', 'figure');
end

o = struct('settingsFile', '', 'ensembleFile', '', 'outDir', '', 'prefix', '');
[S, F] = tdgp3d_resolveRun(runDir, o);

if (ischar(sel) || isstring(sel)) && strcmpi(sel, 'best')
    sel = local_best(runDir, S, F);                        % then read like an iteration number
end
if ischar(sel) || isstring(sel)
    if ~strcmpi(sel, 'last')
        error('tdgp3d_plotNodes:sel', 'sel must be ''last'', ''best'' or iteration numbers.');
    end
    E = tdgp3d_readEnsemble(F.ensemble, 'lastOnly', true);
elseif isnumeric(sel) && ~isempty(sel)
    E = tdgp3d_readEnsemble(F.ensemble, 'burnin', min(sel(:)) - 1);
    keep = ismember(E.imcmc, sel(:));
    E.imcmc = E.imcmc(keep);
    E.nodes = E.nodes(keep);
else
    error('tdgp3d_plotNodes:sel', 'sel must be ''last'', ''best'' or iteration numbers.');
end
if isempty(E.imcmc)
    error('tdgp3d_plotNodes:noSamples', 'No matching samples in %s', F.ensemble);
end
N    = vertcat(E.nodes{:});                                % [x y z log10(rho)], m
L    = [S.lambda_x S.lambda_y S.lambda_z] + S.lambda_increase * max(N(:,3), 0);   % length scales at the nodes
one  = isscalar(E.imcmc);

delete(findobj(figs, 'tag', 'tdgpNodes'));
h = gobjects(0, 1);
nShown = zeros(0, 1);                                      % nodes shown in each panel
for ax = findobj(figs, 'Type', 'axes').'
    switch get(ax, 'Tag')
        case 'tdgp3dMap'                                   % map at depth d: east y across, north x up
            d  = get(ax, 'UserData');
            in = abs(N(:,3) - d) <= L(:,3);
            hx = N(in,2);  hy = N(in,1);
        case 'tdgp3dSectionX'                              % section along x at y = y0
            in = abs(N(:,2) - get(ax, 'UserData')) <= L(:,2);
            hx = N(in,1);  hy = N(in,3);
        case 'tdgp3dSectionY'                              % section along y at x = x0
            in = abs(N(:,1) - get(ax, 'UserData')) <= L(:,1);
            hx = N(in,2);  hy = N(in,3);
        otherwise
            continue;
    end
    nShown(end+1) = nnz(in);                               %#ok<AGROW>
    if ~any(in)
        continue;                                          % no node near this panel
    end
    hold(ax, 'on');
    if one
        hl = plot(ax, hx/1e3, hy/1e3, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
                  'MarkerEdgeColor', 'k', 'LineStyle', 'none', 'Tag', 'tdgpNodes');
    else
        hl = plot(ax, hx/1e3, hy/1e3, 'k.', 'MarkerSize', 6, 'LineStyle', 'none', 'Tag', 'tdgpNodes');
    end
    h = [h; hl];                                           %#ok<AGROW>
end
if isempty(nShown)
    error('tdgp3d_plotNodes:fig', 'No panels of tdgp3d_plotModel found in figs: pass the figures it returned.');
end

for f = figs(:).'
    png = get(f, 'UserData');
    if ischar(png) && ~isempty(png)
        [pth, nm] = fileparts(png);
        tdgp3d_savePNG(f, fullfile(pth, [nm '.nodes.png']));
        fprintf('  saved %s\n', fullfile(pth, [nm '.nodes.png']));
    end
end
fprintf('tdgp3d_plotNodes: %d nodes from %d sample(s), iterations %d to %d; %d to %d nodes shown per panel\n', ...
        size(N, 1), numel(E.imcmc), min(E.imcmc), max(E.imcmc), min(nShown), max(nShown));

end


function it = local_best(runDir, S, F)
% Iteration of the beta = 1 chain with the smallest chi-squared misfit, from
% <outputBase>_summary.txt (columns: iteration, chain, logL, chi2, k, acceptance, ...).
f = fullfile(runDir, [S.outputbase '_summary.txt']);
if exist(f, 'file') ~= 2
    error('tdgp3d_plotNodes:noSummary', '''best'' needs the summary file %s', f);
end
Sm = readmatrix(f, 'FileType', 'text');
[chi2, i] = min(Sm(:,4));
it  = Sm(i,1);
msg = sprintf('best-fitting sample: iteration %d, chi2 %.4g', it, chi2);
if ~isempty(F.data)
    D   = tdgp3d_readModEM(F.data);
    msg = sprintf('%s, RMS %.3f', msg, sqrt(chi2 / (2*numel(D.period))));
end
fprintf('tdgp3d_plotNodes: %s, k = %d\n', msg, Sm(i,5));
end
