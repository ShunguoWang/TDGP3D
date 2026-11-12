function R = tdgp3d_readSites(file)
% TDGP3D_READSITES  MT site codes and positions from a ModEM data file.
%
%   R = tdgp3d_readSites('data.dat')
%
% Every data row is  Period Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error;
% each site is listed once, in the order it first appears.
%
% Output struct: R.code {n x 1}, R.x [n x 1] (north, m), R.y [n x 1] (east, m).

rows = regexp(fileread(file), '\r\n|\n|\r', 'split');
code = {};
x    = zeros(0, 1);
y    = zeros(0, 1);
for i = 1:numel(rows)
    L = strtrim(rows{i});
    if isempty(L) || L(1) == '#' || L(1) == '>'
        continue;
    end
    t = regexp(L, '\s+', 'split');
    if numel(t) < 11 || any(strcmp(code, t{2}))
        continue;
    end
    code{end+1, 1} = t{2};                   %#ok<AGROW>
    x(end+1, 1)    = str2double(t{5});       %#ok<AGROW>
    y(end+1, 1)    = str2double(t{6});       %#ok<AGROW>
end

R.code = code;
R.x    = x;
R.y    = y;

end
