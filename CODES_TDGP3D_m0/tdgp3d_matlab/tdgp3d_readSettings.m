function s = tdgp3d_readSettings(file, base)
% TDGP3D_READSETTINGS  Read a 3D TDGP tdgp.settings Fortran namelist into a struct.
%
%   s = tdgp3d_readSettings('tdgp.settings')
%   s = tdgp3d_readSettings('tdgp3d_out_domain.txt', s)   % overlay a second file onto s
%
% Field names are the namelist variable names in lower case, e.g. s.lambda_x,
% s.lambda_y, s.lambda_z, s.datafile, s.outputbase. Anything the file leaves out
% keeps the built-in default of tdgp3d_com.f90 (or, when a struct 'base' is
% given, the value it already has), which is also what the Fortran code does.
%
% The domain for the node locations (xmin ... zmax) is derived at run time; the
% domain a run used is in <outputBase>_domain.txt, which tdgp3d_resolveRun overlays.
% sigma_m is fixed in tdgp3d_com.f90.

if nargin >= 2 && isstruct(base)
    s = base;
else
    % Defaults (copied from tdgp3d_com.f90 -- keep in sync if that file changes):
    s = struct( ...
        'kmin', 8, 'kmax', 150, 'logrho_min', -1, 'logrho_max', 5, 'rho_ref', 100, ...
        'xy_margin_frac', 0.5, 'xy_margin_min', 5000, 'z_depth_min', 0, ...
        'lambda_x', 1000, 'lambda_y', 1000, 'lambda_z', 500, 'lambda_increase', 0, ...
        'pertsdsc', 30, 'dtlog', 1.15, 'n_chains', 8, 'nproc_per_team', 6, ...
        'nchain', 100000, 'datafile', 'data.dat', 'outputbase', 'tdgp3d_out', ...
        'sigma_m', 0.05, ...
        'xmin', -3000, 'xmax', 3000, 'ymin', -3000, 'ymax', 3000, 'zmin', 0, 'zmax', 3000);
end

if exist(file, 'file') ~= 2
    error('tdgp3d_readSettings:noFile', 'Settings file not found: %s', file);
end

rows = regexp(fileread(file), '\r\n|\n|\r', 'split');
for i = 1:numel(rows)
    L  = rows{i};
    ic = find(L == '!', 1);              % strip Fortran comments
    if ~isempty(ic)
        L = L(1:ic-1);
    end
    tok = regexp(L, '([A-Za-z]\w*)\s*=\s*(''[^'']*''|"[^"]*"|[^,\s/]+)', 'tokens');
    for j = 1:numel(tok)
        key = lower(tok{j}{1});
        val = strtrim(tok{j}{2});
        if val(1) == '''' || val(1) == '"'
            s.(key) = val(2:end-1);                          % quoted string
        else
            num = str2double(regexprep(val, '[dD]', 'e'));   % Fortran 1.0d3 -> 1.0e3
            if isnan(num)
                s.(key) = val;
            else
                s.(key) = num;
            end
        end
    end
end

end
