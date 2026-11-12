function D = tdgp3d_readModEM(file)
% TDGP3D_READMODEM  Data rows of a ModEM data file, as TDGP3D reads them (tdgp3d_modem.f90).
%
%   D = tdgp3d_readModEM('data.dat')
%
% Row r of the file gives data 2r-1 (real part) and 2r (imaginary part) of the TDGP3D data
% vector -- the order of the predicted data in <outputBase>_responses.txt. Supported blocks:
% Full_Impedance, Off_Diagonal_Impedance and Full_Vertical_Components.
%
% Output struct (n = number of data rows):
%   D.period [n x 1] (s), D.code {n x 1}, D.x [n x 1] (north, m), D.y [n x 1] (east, m)
%   D.comp   {n x 1}  'ZXX' 'ZXY' 'ZYX' 'ZYY' 'TX' 'TY'
%   D.re, D.im, D.err [n x 1]  as in the file
%   D.toOhm  [n x 1]  impedance in Ohm = value * toOhm (1 for the tipper)
%   D.sgn    [n x 1]  +1 for exp(+i omega t) rows, -1 for exp(-i omega t) rows

mu0  = 4e-7 * pi;
rows = regexp(fileread(file), '\r\n|\n|\r', 'split');
n    = 0;
for i = 1:numel(rows)
    L = strtrim(rows{i});
    if ~isempty(L) && L(1) ~= '#' && L(1) ~= '>'
        n = n + 1;
    end
end
D.period = zeros(n,1);  D.code = cell(n,1);  D.x = zeros(n,1);  D.y = zeros(n,1);
D.comp   = cell(n,1);   D.re   = zeros(n,1); D.im = zeros(n,1); D.err = zeros(n,1);
D.toOhm  = ones(n,1);   D.sgn  = ones(n,1);

r = 0; nhead = 0; inData = false; btype = ''; sgn = 1; toOhm = 1;
for i = 1:numel(rows)
    L = strtrim(rows{i});
    if isempty(L)
        continue;
    end
    if L(1) == '#'
        nhead = 0; inData = false;
        continue;
    end
    if L(1) == '>'
        if inData
            nhead = 0; inData = false;
        end
        nhead = nhead + 1;
        s = lower(strtrim(L(2:end)));
        switch nhead
            case 1
                if     startsWith(s, 'full_impedance'),           btype = 'Z';
                elseif startsWith(s, 'off_diagonal_impedance'),   btype = 'Z';
                elseif startsWith(s, 'full_vertical_components'), btype = 'T';
                else
                    error('tdgp3d_readModEM:type', '%s: unsupported data type "%s"', file, strtrim(L(2:end)));
                end
            case 2
                if contains(s, '+'), sgn = 1; elseif contains(s, '-'), sgn = -1;
                else, error('tdgp3d_readModEM:sign', '%s: cannot read the sign convention "%s"', file, L);
                end
            case 3
                if btype == 'T',                                   toOhm = 1;
                elseif contains(s, '[mv/km]/[nt]'),                toOhm = mu0 * 1000;
                elseif contains(s, '[v/m]/[t]'),                   toOhm = mu0;
                elseif contains(s, '[v/m]/[a/m]') || contains(s, 'ohm'), toOhm = 1;
                else, error('tdgp3d_readModEM:units', '%s: unknown impedance units "%s"', file, L);
                end
        end
        continue;
    end
    inData = true;
    t = regexp(L, '\s+', 'split');
    if numel(t) < 11
        error('tdgp3d_readModEM:row', '%s: cannot read the data row: %s', file, L);
    end
    r = r + 1;
    c = upper(t{8});
    if strcmp(c, 'TZX'), c = 'TX'; end
    if strcmp(c, 'TZY'), c = 'TY'; end
    D.period(r) = str2double(t{1});
    D.code{r}   = t{2};
    D.x(r)      = str2double(t{5});
    D.y(r)      = str2double(t{6});
    D.comp{r}   = c;
    D.re(r)     = str2double(t{9});
    D.im(r)     = str2double(t{10});
    D.err(r)    = str2double(t{11});
    D.toOhm(r)  = toOhm;
    D.sgn(r)    = sgn;
end

end
