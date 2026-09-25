function RHS = RHSbracketmatrix(nv_omega, nv_theta, nf_omega, nf_theta, V, F)

% ============================================================
% Cached RHS bracket matrix
% - Derivative operators are built ONCE
% - Subsequent calls reuse them
% ============================================================

persistent Dcache keycache

% ------------------------------------------------------------
% Build cache key
% ------------------------------------------------------------
key = {nv_omega, nv_theta, nf_omega, nf_theta};

if isempty(Dcache) || ~isequal(key, keycache)

    c = numel(nf_omega);
    d = numel(nf_theta);

    % ---------- omega operators ----------
if c~=0
    Dwlf = cell(c,1);
    Dwlv = cell(c,1);
    parfor l = 1:c
        Dwlf{l} = D_omega(l, nf_omega, nf_theta);
        Dwlv{l} = D_omega(l, nv_omega, nv_theta);
    end
end
% ---------- theta operators ----------
Dtlf = cell(d,1);
Dtlv = cell(d,1);

parfor l = 1:d
    Dtlf{l} = D_theta(l, nf_omega, nf_theta);
    Dtlv{l} = D_theta(l, nv_omega, nv_theta);
end

% ---------- pack into struct ----------
if c~=0
    Dcache.Dwlf = Dwlf;
    Dcache.Dwlv = Dwlv;
end
Dcache.Dtlf = Dtlf;
Dcache.Dtlv = Dtlv;

keycache = key;
end

% ------------------------------------------------------------
% Assemble RHS
% ------------------------------------------------------------
c = numel(nf_omega);
d = numel(nf_theta);

rsize = prod(nv_theta+1) * prod(nv_omega+1) * ...
        prod(nf_theta+1) * prod(nf_omega+1);

RHS = sparse(rsize, rsize);

if c~=0
    for l = 1:c
        RHS = RHS + kron(V, Dcache.Dwlf{l}' * F{l} + F{l} * Dcache.Dwlf{l}) ...
              - kron(Dcache.Dwlv{l}' * V + V * Dcache.Dwlv{l}, F{l});
    end
end

for l = 1:d
    RHS = RHS + kron(V, Dcache.Dtlf{l}' * F{c+l} + F{c+l} * Dcache.Dtlf{l}) ...
              - kron(Dcache.Dtlv{l}' * V + V * Dcache.Dtlv{l}, F{c+l});
end

end
