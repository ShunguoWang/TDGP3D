function E = tdgp3d_readEnsemble(file, varargin)
% TDGP3D_READENSEMBLE  Read samples from a 3D TDGP ensemble file (<outputBase>_ensemble.txt).
%
%   E = tdgp3d_readEnsemble(file)
%   E = tdgp3d_readEnsemble(file, 'burnin',2000, 'thin',10, 'maxSamples',1000)
%   E = tdgp3d_readEnsemble(file, 'lastOnly',true)
%
% Each row of the file is one sample (one model from the chain at beta = 1):
%       imcmc  chain_id  k   x1 y1 z1 logrho1   x2 y2 z2 logrho2  ...  xk yk zk logrhok  [m0]
% with x north, y east, z depth (m) and logrho = log10(ohm-m). k, the number of GP
% nodes, changes from row to row. The nodes are the control points of a smooth
% field; tdgp3d_gpMean turns them into resistivity on the forward mesh. m0, the
% sample's background (GP prior mean and forward half-space), ends the row only in
% runs with sample_m0; otherwise every sample has the background of MT3D_V6.dat.
%
% Options (name/value):
%   'burnin'      discard rows with iteration (imcmc) <= burnin          [0]
%   'thin'        keep every thin-th row that passes the burn-in         [1]
%   'maxSamples'  stop after this many kept rows (memory guard)          [Inf]
%   'lastOnly'    keep only the last row                                 false
%
% Output struct E (n = number of rows kept):
%   E.imcmc [n x 1], E.chain [n x 1], E.k [n x 1]
%   E.nodes {n x 1}  each cell is a [k x 4] matrix: [x y z log10(rho)]
%   E.m0    [n x 1]  the sample's background log10(rho) (sample_m0 runs); NaN where the
%                    row has none (the background of MT3D_V6.dat applies)
%   E.nBad           rows skipped because they could not be parsed (e.g. a line cut short)
%   E.hitMaxSamples  true if reading stopped early because of 'maxSamples'

p = inputParser;
p.FunctionName = 'tdgp3d_readEnsemble';
addParameter(p, 'burnin',     0,     @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'thin',       1,     @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'maxSamples', Inf,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'lastOnly',   false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
o = p.Results;

fid = fopen(file, 'r');
if fid < 0
    error('tdgp3d_readEnsemble:open', 'Cannot open ensemble file: %s', file);
end
closer = onCleanup(@() fclose(fid));

% A file that does not end with a newline had its last row cut short (e.g. the
% job was killed while writing); that row must not be used.
endsWithNewline = true;
if fseek(fid, -1, 'eof') == 0
    endsWithNewline = isequal(fread(fid, 1, '*char'), newline);
end
fseek(fid, 0, 'bof');

cap   = 1024;
imc   = zeros(cap, 1);
chn   = zeros(cap, 1);
kk    = zeros(cap, 1);
nod   = cell(cap, 1);
m0v   = nan(cap, 1);
n     = 0;          % rows kept
nElig = 0;          % rows past the burn-in
nBad  = 0;          % rows that could not be parsed
hitCap   = false;
lastLine = '';
prevLine = '';

while true
    L = fgetl(fid);
    if ~ischar(L)
        break;
    end
    if ~endsWithNewline && feof(fid)
        if ~isempty(strtrim(L))
            nBad = nBad + 1;
        end
        break;
    end
    h3 = sscanf(L, '%f', 3);                     % imcmc, chain_id, k
    if numel(h3) < 3
        if ~isempty(strtrim(L))
            nBad = nBad + 1;
        end
        continue;
    end
    if h3(1) <= o.burnin
        continue;
    end
    if o.lastOnly
        prevLine = lastLine;
        lastLine = L;
        continue;
    end
    nElig = nElig + 1;
    if mod(nElig - 1, o.thin) ~= 0
        continue;
    end
    v  = sscanf(L, '%f');
    nv = 3 + 4*h3(3);                            % values of a row without m0
    if numel(v) ~= nv && numel(v) ~= nv + 1
        nBad = nBad + 1;
        continue;
    end
    n = n + 1;
    if n > cap
        cap = 2*cap;
        imc(cap) = 0;
        chn(cap) = 0;
        kk(cap)  = 0;
        nod{cap} = [];
        m0v(cap) = NaN;
    end
    imc(n) = v(1);
    chn(n) = v(2);
    kk(n)  = v(3);
    nod{n} = reshape(v(4:nv), 4, v(3)).';        % [k x 4] = [x y z logrho]
    m0v(n) = local_m0(v, nv);
    if n >= o.maxSamples
        hitCap = true;
        break;
    end
end

if o.lastOnly
    cand = {lastLine, prevLine};                  % the very last line may be cut short
    for c = 1:2
        if isempty(cand{c})
            continue;
        end
        v = sscanf(cand{c}, '%f');
        if numel(v) >= 3 && (numel(v) == 3 + 4*v(3) || numel(v) == 4 + 4*v(3))
            nv = 3 + 4*v(3);
            n = 1;
            imc(1) = v(1);
            chn(1) = v(2);
            kk(1)  = v(3);
            nod{1} = reshape(v(4:nv), 4, v(3)).';
            m0v(1) = local_m0(v, nv);
            break;
        else
            nBad = nBad + 1;
        end
    end
end

if hitCap
    warning('tdgp3d_readEnsemble:maxSamples', ...
        ['Stopped after maxSamples = %d kept rows (at iteration %d) but the file has more. ' ...
         'Increase ''thin'' (or ''maxSamples'') to cover the whole run.'], o.maxSamples, imc(n));
end
if nBad > 0
    warning('tdgp3d_readEnsemble:badRows', ...
        '%d row(s) could not be parsed and were skipped (a line cut short?).', nBad);
end

E.imcmc = imc(1:n);
E.chain = chn(1:n);
E.k     = kk(1:n);
E.nodes = nod(1:n);
E.m0    = m0v(1:n);
E.nBad  = nBad;
E.hitMaxSamples = hitCap;

end

function m0 = local_m0(v, nv)
% The background m0 ending a row of a sample_m0 run, NaN for a row without it.
if numel(v) == nv + 1
    m0 = v(end);
else
    m0 = NaN;
end
end
