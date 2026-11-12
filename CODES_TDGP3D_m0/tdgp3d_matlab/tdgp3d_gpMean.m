function [mu, ok, alpha] = tdgp3d_gpMean(nodes, xs, ys, zs, lambda, sigma_m, mean0, inc)
% TDGP3D_GPMEAN  GP conditional mean of log10 resistivity at given points (3D).
%
%   [mu, ok] = tdgp3d_gpMean(nodes, xs, ys, zs, [lambda_x lambda_y lambda_z], sigma_m, logrho_bg)
%   [~, ok, alpha] = tdgp3d_gpMean(nodes, [], [], [], lambda, sigma_m, logrho_bg)   % weights only
%   [mu, ok] = tdgp3d_gpMean(..., logrho_bg, lambda_increase)    % length scales grow with depth
%
%   nodes       [k x 4]  columns x, y, z, log10(rho) of the k GP nodes
%   xs, ys, zs  [n x 1]  points to evaluate (e.g. the earth nodes of the forward mesh)
%   lambda      [1 x 3]  GP length scales lambda_x, lambda_y, lambda_z (m), from tdgp.settings
%   sigma_m              nugget (0.05, fixed in tdgp3d_com.f90)
%   mean0                GP prior mean = log10 background resistivity of MT3D_V6.dat
%   inc                  (optional) lambda_increase from tdgp.settings [m per m of depth]:
%                        every length scale is lambda + inc*depth at each node and point
%                        (depth = max(z, 0)); 0 or left out = constant length scales
%
%   mu     [n x 1]  log10(rho) at the points;  ok = false if the node kernel matrix is
%                   not positive definite (mu = NaN)
%   alpha  [k x 1]  GP weights Km \ (g - mean0), the same for any points:
%                   mu = mean0 + K(x,x*)' * alpha  ([] if ok is false)
%
% Same formulas as GP_PREDICT_MEAN3D in tdgp3d_gp.f90:
%   K(a,b) = exp( -0.5*[ ((xa-xb)/lambda_x)^2 + ((ya-yb)/lambda_y)^2 + ((za-zb)/lambda_z)^2 ] )
%   Km     = K(x,x) + sigma_m^2 * I
%   mu     = mean0 + K(x,x*)' * ( Km \ (g - mean0) )
% With inc > 0 the kernel is that of GP_KERNEL3D_INC, with l_c = lambda_c + inc*max(z, 0):
%   K(a,b) = PROD_c sqrt(2 l_c(a) l_c(b)/(l_c(a)^2 + l_c(b)^2)) * exp(-SUM_c (c_a-c_b)^2/(l_c(a)^2 + l_c(b)^2))

if nargin >= 8 && inc > 0                               % length scales grow with depth
    [mu, ok, alpha] = local_gpMeanInc(nodes, [xs(:) ys(:) zs(:)], lambda, sigma_m, mean0, inc);
    return;
end
xs = xs(:);
ys = ys(:);
zs = zs(:);
n  = numel(xs);
k  = size(nodes, 1);

if k == 0
    mu = mean0 * ones(n, 1);        % no nodes: the prior mean everywhere
    ok = true;
    alpha = zeros(0, 1);
    return;
end

x = nodes(:,1) / lambda(1);
y = nodes(:,2) / lambda(2);
z = nodes(:,3) / lambda(3);
g = nodes(:,4) - mean0;

Km = exp(-0.5 * ((x - x.').^2 + (y - y.').^2 + (z - z.').^2)) + sigma_m^2 * eye(k);
[R, flag] = chol(Km);                                   % Km = R'*R
if flag ~= 0
    mu = nan(n, 1);
    ok = false;
    alpha = [];
    return;
end
alpha = R \ (R' \ g);

mu    = zeros(n, 1);
chunk = max(1, floor(2e6 / k));                         % bounds the k x chunk kernel block
xs = xs / lambda(1);
ys = ys / lambda(2);
zs = zs / lambda(3);
for i0 = 1:chunk:n
    i1 = min(n, i0 + chunk - 1);
    Ks = exp(-0.5 * ((x - xs(i0:i1).').^2 + (y - ys(i0:i1).').^2 + (z - zs(i0:i1).').^2));
    mu(i0:i1) = Ks.' * alpha;
end
mu = mean0 + mu;
ok = true;

end


function [mu, ok, alpha] = local_gpMeanInc(nodes, pts, lambda, sigma_m, mean0, inc)
% tdgp3d_gpMean with length scales lambda + inc*depth; pts [n x 3] x, y, z (m).
n = size(pts, 1);
k = size(nodes, 1);
if k == 0
    mu = mean0 * ones(n, 1);
    ok = true;
    alpha = zeros(0, 1);
    return;
end
[R, flag] = chol(local_kernelInc(nodes(:,1:3), nodes(:,1:3), lambda, inc) + sigma_m^2 * eye(k));
if flag ~= 0
    mu = nan(n, 1);
    ok = false;
    alpha = [];
    return;
end
alpha = R \ (R' \ (nodes(:,4) - mean0));
mu    = zeros(n, 1);
chunk = max(1, floor(2e6 / k));
for i0 = 1:chunk:n
    i1 = min(n, i0 + chunk - 1);
    mu(i0:i1) = local_kernelInc(nodes(:,1:3), pts(i0:i1,:), lambda, inc).' * alpha;
end
mu = mean0 + mu;
ok = true;
end


function K = local_kernelInc(p1, p2, lambda, inc)
% K(i,j) between the points p1(i,:) and p2(j,:), [x y z] in m (GP_KERNEL3D_INC).
l1 = lambda(:).' + inc * max(p1(:,3), 0);              % length scales at the points
l2 = lambda(:).' + inc * max(p2(:,3), 0);
K  = 1;
E  = 0;
for c = 1:3
    s = l1(:,c).^2 + l2(:,c).'.^2;
    K = K .* sqrt(2 * l1(:,c) .* l2(:,c).' ./ s);
    E = E - (p1(:,c) - p2(:,c).').^2 ./ s;
end
K = K .* exp(E);
end
