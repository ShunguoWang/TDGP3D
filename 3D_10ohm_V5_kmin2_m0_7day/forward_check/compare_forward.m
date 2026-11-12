% COMPARE_FORWARD  Forward responses of the test model from MTSphere and from V6 (the inversion's forward
% code) at the 361 stations of 3D_10ohm_V3, at 10, 1 and 0.1 Hz.
%
%   Run from this folder:  compare_forward
%
% Model: 10 ohm-m sphere, radius 500 m, centre (0, 0, 1000 m), in a 1000 ohm-m half-space.
%   MTSphere (semi-analytic):      ../data_noisefree.dat
%   V6 on the inversion mesh:      ../true_model_V6_response.txt (the sphere given at the nodes of MT3D_V6.dat, as the
%                                  inversion passes its models; same rows and units as data_noisefree.dat)
%   errors used by the inversion:  ../data.dat (5% of sqrt|Zxy Zyx| for Z, 0.01 for the tipper)
% Figure 1: apparent resistivity and phase of Zxy and Zyx along the profiles y = 0 and x = 0.
% Figure 2: maps of V6 - MTSphere: apparent resistivity (% of MTSphere) and phase (degrees).
% Prints a table per frequency. Phases for exp(+i omega t): Zxy near 45 deg, Zyx shown as the phase of -Zyx.
% Writes compare_forward_profiles.png/.pdf and compare_forward_maps.png/.pdf to outDir (default: this folder).

addpath('/Users/shunguow/Desktop/MT3D_Bayesian/CODES_TDGP3D/tdgp3d_matlab');
here = fileparts(mfilename('fullpath'));
if ~exist('outDir', 'var'), outDir = here; end
runDir = fullfile(here, '..');
R = 500;  mu0 = 4e-7*pi;

D = tdgp3d_readModEM(fullfile(runDir, 'data_noisefree.dat'));
E = tdgp3d_readModEM(fullfile(runDir, 'data.dat'));
L = regexp(fileread(fullfile(runDir, 'true_model_V6_response.txt')), '\r\n|\n|\r', 'split');
L = L(~cellfun(@isempty, L) & ~startsWith(L, '%'));
p = str2double(L(:));
n = numel(D.re);
if numel(p) ~= 2*n || numel(E.re) ~= n, error('the files do not have the same rows'); end
zm = (D.re + 1i*D.im) .* D.toOhm;              zm(D.sgn < 0) = conj(zm(D.sgn < 0));   % Ohm, exp(+i omega t)
zv = (p(1:2:end) + 1i*p(2:2:end)) .* D.toOhm;  zv(D.sgn < 0) = conj(zv(D.sgn < 0));
er = E.err .* E.toOhm;
f  = 1 ./ D.period;
fr = [10 1 0.1];
kf = zeros(n, 1);
for k = 1:3, kf(abs(f - fr(k)) <= 1e-6*fr(k)) = k; end
isXY = strcmp(D.comp, 'ZXY');  isYX = strcmp(D.comp, 'ZYX');  isZ = startsWith(D.comp, 'Z');
rhoa = @(z, fq) abs(z).^2 ./ (2*pi*fq*mu0);
phs  = @(z) atan2d(imag(z), real(z));
wrap = @(a) mod(a + 180, 360) - 180;

% --- table ---------------------------------------------------------------------------------------------
fprintf('\nV6 - MTSphere at %d stations (errors of the inversion data)\n', numel(unique(D.code)));
fprintf('  f (Hz) | rho_a xy, yx: RMS / max (%%) | phase xy, yx: RMS / max (deg) | RMS in errors: Z, T, all\n');
for k = 1:3
    dR = []; dP = [];
    for c = {isXY, isYX}
        r = c{1} & kf == k;  s = 1 - 2*isYX(find(r, 1));
        dR = [dR 100*(rhoa(zv(r), f(r)) ./ rhoa(zm(r), f(r)) - 1)];  %#ok<AGROW>
        dP = [dP wrap(phs(s*zv(r)) - phs(s*zm(r)))];                   %#ok<AGROW>
    end
    rz = isZ & kf == k;  rt = ~isZ & kf == k;
    nz = [real(zv(rz) - zm(rz)); imag(zv(rz) - zm(rz))] ./ [er(rz); er(rz)];
    nt = [real(zv(rt) - zm(rt)); imag(zv(rt) - zm(rt))] ./ [er(rt); er(rt)];
    fprintf('  %6g |  %4.2f %4.2f / %4.1f %4.1f     |  %5.2f %5.2f / %5.2f %5.2f         |  %5.3f %5.3f %5.3f\n', fr(k), ...
            sqrt(mean(dR.^2)), max(abs(dR)), sqrt(mean(dP.^2)), max(abs(dP)), ...
            sqrt(mean(nz.^2)), sqrt(mean(nt.^2)), sqrt(mean([nz; nt].^2)));
end
nz = [real(zv - zm); imag(zv - zm)] ./ [er; er];
fprintf('  all    | RMS in the errors of the inversion: %.3f (the noise added to the data has RMS 1)\n\n', sqrt(mean(nz.^2)));

% --- figure 1: profiles y = 0 and x = 0 ------------------------------------------------------------------
cx = [0 0.35 0.75];  cy = [0.85 0.2 0.1];
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 19 16]);
tl  = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
prof = {D.y == 0, D.x == 0};  pos = {D.x, D.y};  plab = {'x, north (km), along y = 0', 'y, east (km), along x = 0'};
for k = 1:3
    for q = 1:2
        for w = 1:2                                         % w = 1: rho_a, w = 2: phase
            ax = nexttile(tl);  hold(ax, 'on');
            for c = 1:2
                r = prof{q} & kf == k & ternary(c == 1, isXY, isYX);
                [xs, o] = sort(pos{q}(r));  ri = find(r);  ri = ri(o);  s = 3 - 2*c;   % s = -1 for Zyx
                if w == 1, ym = rhoa(zm(ri), f(ri)); yv = rhoa(zv(ri), f(ri));
                else,      ym = phs(s*zm(ri));        yv = phs(s*zv(ri)); end
                col = ternary(c == 1, cx, cy);
                plot(ax, xs/1e3, ym, '-', 'Color', col, 'LineWidth', 1.2);
                plot(ax, xs/1e3, yv, 'o', 'Color', col, 'MarkerSize', 3.5, 'LineWidth', 0.8);
            end
            set(ax, 'Box', 'on', 'FontSize', 7, 'XLim', [-1.9 1.9]);
            xline(ax, [-R R]/1e3, ':', 'Color', [0.5 0.5 0.5]);
            if k == 1
                title(ax, sprintf('%s, %s', ternary(w == 1, '\rho_a (\Omega m)', 'phase (\circ)'), ...
                      ternary(q == 1, 'y = 0', 'x = 0')), 'FontSize', 8, 'FontWeight', 'normal');
            end
            if k == 3, xlabel(ax, plab{q}); else, set(ax, 'XTickLabel', []); end
            if q == 1 && w == 1, ylabel(ax, sprintf('%g Hz', fr(k)), 'FontSize', 8, 'FontWeight', 'bold'); end
        end
    end
end
lg = legend(nexttile(tl, 1), {'Z_{xy} MTSphere', 'Z_{xy} V6', 'Z_{yx} MTSphere', 'Z_{yx} V6'}, ...
            'Orientation', 'horizontal', 'FontSize', 7);
lg.Layout.Tile = 'north';
title(tl, 'MTSphere (lines) and V6 on the inversion mesh (circles); dotted: edges of the sphere', 'FontSize', 8);
exportgraphics(fig, fullfile(outDir, 'compare_forward_profiles.png'), 'Resolution', 250);
exportgraphics(fig, fullfile(outDir, 'compare_forward_profiles.pdf'), 'ContentType', 'vector');

% --- figure 2: maps of the differences -----------------------------------------------------------------
ux = unique(D.x);  uy = unique(D.y);
cm = [linspace(0.15, 1, 128)' linspace(0.3, 1, 128)' ones(128, 1); ones(128, 1) linspace(1, 0.25, 128)' linspace(1, 0.15, 128)'];
t  = linspace(0, 2*pi, 200);
M  = cell(4, 3);
for k = 1:3
    for c = 1:2
        sel = kf == k & ternary(c == 1, isXY, isYX);  s = 3 - 2*c;
        A = nan(numel(ux), numel(uy));  B = A;
        for i = find(sel).'
            ix = ux == D.x(i);  iy = uy == D.y(i);
            A(ix, iy) = 100*(rhoa(zv(i), f(i)) / rhoa(zm(i), f(i)) - 1);
            B(ix, iy) = wrap(phs(s*zv(i)) - phs(s*zm(i)));
        end
        M{c, k} = A;  M{c + 2, k} = B;
    end
end
lim = [max(abs([M{1:2, :}]), [], 'all') max(abs([M{3:4, :}]), [], 'all')];
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 21]);
tl  = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rl  = {'\Delta\rho_{a,xy} (%)', '\Delta\rho_{a,yx} (%)', '\Delta\phi_{xy} (\circ)', '\Delta\phi_{yx} (\circ)'};
for row = 1:4
    for k = 1:3
        ax = nexttile(tl);  hold(ax, 'on');
        imagesc(ax, uy/1e3, ux/1e3, M{row, k});
        plot(ax, R*cos(t)/1e3, R*sin(t)/1e3, 'k--', 'LineWidth', 0.8);
        L2 = lim(1 + (row > 2));
        set(ax, 'YDir', 'normal', 'Box', 'on', 'FontSize', 7, 'DataAspectRatio', [1 1 1], 'CLim', [-L2 L2], ...
            'XLim', [-1.9 1.9], 'YLim', [-1.9 1.9], 'Layer', 'top');
        colormap(ax, cm);
        if row == 1, title(ax, sprintf('%g Hz', fr(k)), 'FontSize', 8, 'FontWeight', 'normal'); end
        if row == 4, xlabel(ax, 'y, east (km)'); else, set(ax, 'XTickLabel', []); end
        if k == 1, ylabel(ax, {rl{row}, 'x, north (km)'}); else, set(ax, 'YTickLabel', []); end
        if k == 3, cb = colorbar(ax); cb.FontSize = 7; end
    end
end
title(tl, 'V6 - MTSphere at the 361 stations; dashed: outline of the sphere', 'FontSize', 8);
exportgraphics(fig, fullfile(outDir, 'compare_forward_maps.png'), 'Resolution', 250);
exportgraphics(fig, fullfile(outDir, 'compare_forward_maps.pdf'), 'ContentType', 'vector');

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
