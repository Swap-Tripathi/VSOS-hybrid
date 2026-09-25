function RHS = RHSbracketmatrix_LaSalle_local(nv_omega, nv_theta, nf_omega, nf_theta, V, F, Q)

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
Dwlv = cell(c,1);

parfor l = 1:c
    Dwlv{l} = D_omega(l, nv_omega, nv_theta);
end

% ---------- theta operators ----------
Dtlv = cell(d,1);

parfor l = 1:d
    Dtlv{l} = D_theta(l, nv_omega, nv_theta);
end

% ---------- pack into struct ----------
Dcache.Dwlv = Dwlv;
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

for l = 1:c
    RHS = RHS - kron(Dcache.Dwlv{l}' * V + V * Dcache.Dwlv{l}, F{l});
end

for l = 1:d
    RHS = RHS - kron(Dcache.Dtlv{l}' * V + V * Dcache.Dtlv{l}, F{c+l});
end

end
