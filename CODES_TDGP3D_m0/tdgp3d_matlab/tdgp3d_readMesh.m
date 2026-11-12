function M = tdgp3d_readMesh(file)
% TDGP3D_READMESH  Node coordinates of the V6 forward mesh from MT3D_V6.dat.
%
%   M = tdgp3d_readMesh('MT3D_V6.dat')
%
% Reads the file the way V6's READ_PARA does (one value per line after the
% label; the dz list may span several lines and may use the r*c repeat form)
% and builds the node coordinates the way V6's MODEL does:
%   x(i) = x0 + (i-1)*dx,  y(j) = y0 + (j-1)*dy
%   z: node jf (number_grid_in_air) at z = 0; every element has a middle node and
%      is twice its dz value thick; air elements go up from jf, earth elements down.
%
% Output struct:
%   M.x [nx x 1], M.y [ny x 1], M.z [nz x 1]   node coordinates (m), z positive down
%   M.jf          index of the surface node (z = 0); earth nodes are M.z(M.jf:end)
%   M.dx, M.dy, M.dzin [ke x 1] (the dz list as written), M.sigma0 (background, S/m),
%   M.logrho_bg = log10(1/sigma0), M.itermax

rows = regexp(fileread(file), '\r\n|\n|\r', 'split');
rows = rows(~cellfun(@(r) isempty(strtrim(r)), rows));   % list-directed reads skip empty lines
r = 1;
r = r + 2;                                   % title, separator
r = r + 3;                                   % freqmin, freqmax, nfreq (not used)
r = r + 1;                                   % separator
nx = local_value(rows, r);   r = r + 1;
ny = local_value(rows, r);   r = r + 1;
nz = local_value(rows, r);   r = r + 1;
dx = local_value(rows, r);   r = r + 1;
dy = local_value(rows, r);   r = r + 1;
jf = local_value(rows, r);   r = r + 1;
r = r + 2;                                   % separator, dz comment line
ke = (nz - 1) / 2;
dzin = [];
while numel(dzin) < ke
    if r > numel(rows)
        error('tdgp3d_readMesh:format', '%s: the dz list has fewer than (nnz-1)/2 = %d values', file, ke);
    end
    dzin = [dzin; local_numbers(rows{r})];   %#ok<AGROW>
    r = r + 1;
end
dzin = dzin(1:ke);
r = r + 1;                                   % separator
x0 = local_value(rows, r);   r = r + 1;
y0 = local_value(rows, r);   r = r + 1;
r = r + 1;                                   % separator
sigma0  = local_value(rows, r);   r = r + 1;
r = r + 1;                                   % anomal_conductivity (not used)
itermax = local_value(rows, r);

if mod(jf - 1, 2) ~= 0 || jf < 1 || jf > nz
    error('tdgp3d_readMesh:format', '%s: number_grid_in_air = %g must be odd and at most nnz', file, jf);
end
nelair = (jf - 1) / 2;
z = zeros(nz, 1);
z(jf) = 0;
for i = nelair:-1:1                          % air elements, from the surface up
    d = 2 * dzin(i);
    z(2*i)   = z(2*i+1) - d/2;
    z(2*i-1) = z(2*i+1) - d;
end
for i = nelair+1:ke                          % earth elements, from the surface down
    d = 2 * dzin(i);
    z(2*i)   = z(2*i-1) + d/2;
    z(2*i+1) = z(2*i-1) + d;
end

M.x         = x0 + (0:nx-1).' * dx;
M.y         = y0 + (0:ny-1).' * dy;
M.z         = z;
M.jf        = jf;
M.dx        = dx;
M.dy        = dy;
M.dzin      = dzin;
M.sigma0    = sigma0;
M.logrho_bg = log10(1 / sigma0);
M.itermax   = itermax;

end


function v = local_value(rows, r)
% The value on a "label value" line (list-directed: the label is the first item).
t = local_items(rows{r});
if numel(t) < 2
    error('tdgp3d_readMesh:format', 'Expected "label value" on line: %s', rows{r});
end
v = str2double(regexprep(t{2}, '[dD]', 'e'));
if isnan(v)
    error('tdgp3d_readMesh:format', 'Cannot read a number on line: %s', rows{r});
end
end


function v = local_numbers(L)
% All numbers on a line, with Fortran's r*c repeat form (e.g. 75*20.) expanded.
t = local_items(L);
v = zeros(0, 1);
for i = 1:numel(t)
    p = strsplit(t{i}, '*');
    if numel(p) == 2
        v = [v; repmat(str2double(regexprep(p{2}, '[dD]', 'e')), str2double(p{1}), 1)];   %#ok<AGROW>
    else
        v = [v; str2double(regexprep(t{i}, '[dD]', 'e'))];                                  %#ok<AGROW>
    end
end
if any(isnan(v))
    error('tdgp3d_readMesh:format', 'Cannot read the dz values on line: %s', L);
end
end


function t = local_items(L)
t = regexp(strtrim(L), '[\s,]+', 'split');
end
