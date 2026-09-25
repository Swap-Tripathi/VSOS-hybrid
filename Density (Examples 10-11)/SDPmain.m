function [status_main, Vmain, Wmain] = SDPmain( ...
        nv_omega, nv_theta, nf_omega, nf_theta, F, decaypower_p, R)

% ================= DIMENSIONS =================
nw_omega = nv_omega + nf_omega;
nw_theta = nv_theta + nf_theta;
nS0_omega = nv_omega;
nS0_theta = nv_theta;
nS1_omega = nv_omega-2*ones;                  
nS1_theta = nv_theta;
nq1_omega = 2*ones(1,length(nf_omega));
nq1_theta = 0;

Q1=normPowerGramSparse(1,nq1_omega,nq1_theta);
Q1(1,1) = Q1(1,1) - R^2;

c = numel(nw_omega);
d = numel(nw_theta);


% ================= INDEX SETS =================

% ---- Theta indices k ----
grids = cell(1,d);
for j = 1:d
    grids{j} = -nw_theta(j) : nw_theta(j);
end
[Kt{1:d}] = ndgrid(grids{:});

STheta = zeros(numel(Kt{1}), d);
for j = 1:d
    STheta(:,j) = Kt{j}(:);
end
STheta = STheta.';

% ---- Omega indices eta ----
if c~=0
    grids = cell(1,c);
    for j = 1:c
        grids{j} = 0 : 2*nw_omega(j);
    end
    [Ko{1:c}] = ndgrid(grids{:});

    SOmega = zeros(numel(Ko{1}), c);
    for j = 1:c
        SOmega(:,j) = Ko{j}(:);
    end
    SOmega = SOmega.';
end

% ---- Theta indices k for V ----
gridsv = cell(1,d);
for j = 1:d
    gridsv{j} = -nv_theta(j) : nv_theta(j);
end
[Ktv{1:d}] = ndgrid(gridsv{:});

SThetav = zeros(numel(Ktv{1}), d);
for j = 1:d
    SThetav(:,j) = Ktv{j}(:);
end
SThetav = SThetav.';

% ---- Omega indices eta for V----
if c~=0
    gridsv = cell(1,c);
    for j = 1:c
        gridsv{j} = 0 : 2*nv_omega(j);
    end
    [Kov{1:c}] = ndgrid(gridsv{:});

    SOmegav = zeros(numel(Kov{1}), c);
    for j = 1:c
        SOmegav(:,j) = Kov{j}(:);
    end
    SOmegav = SOmegav.';
end

% ================= SELECTION MATRIX =================
X = selection_matrix(nv_omega, nv_theta, nf_omega, nf_theta);

% ================= PRECOMPUTE curlH OPERATORS for W =================
Klist   = STheta(:,1:(end+1)/2).';
if c~=0
    Etalist = SOmega.';
    nEta = size(Etalist,1);
end
nK   = size(Klist,1);
if c==0
    nConstr = nK;
else
    nConstr = nK * nEta;
end
ConstrData = cell(nConstr,1);

if c~=0
    parfor idx = 1:nConstr
        [iK,iEta] = ind2sub([nK,nEta],idx);
        k   = Klist(iK,:).';
        eta = Etalist(iEta,:).';
        ConstrData{idx} = curlyH_operator(nw_omega, nw_theta, eta, k);
    end
else
    parfor idx = 1:nConstr
        [iK,iEta] = ind2sub([nK,1],idx);
        k   = Klist(iK,:).';
        ConstrData{idx} = curlyH_operator(nw_omega, nw_theta, [], k);
    end
end

% ================= PRECOMPUTE curlH OPERATORS for V =================
Klistv   = SThetav(:,1:(end+1)/2).';
if c~=0
    Etalistv = SOmegav.';
    nEtav = size(Etalistv,1);
end
nKv   = size(Klistv,1);
if c==0
    nConstrv = nKv;
else
    nConstrv = nKv * nEtav;
end
ConstrDatav = cell(nConstrv,1);

if c~=0
    parfor idx = 1:nConstrv
        [iKv,iEtav] = ind2sub([nKv,nEtav],idx);
        kv   = Klistv(iKv,:).';
        etav = Etalistv(iEtav,:).';
        ConstrDatav{idx} = curlyH_operator(nv_omega, nv_theta, etav, kv);
    end
else
    parfor idx = 1:nConstrv
        [iKv,iEtav] = ind2sub([nKv,1],idx);
        kv   = Klistv(iKv,:).';
        ConstrDatav{idx} = curlyH_operator(nv_omega, nv_theta, [], kv);
    end
end

% ================= SDP =================
cvx_clear
cvx_solver mosek_2 %MOSEK Version 11.0.29
cvx_precision low
cvx_quiet(false)

cvx_begin sdp

    size_v = prod(nv_theta+1) * prod(nv_omega+1);
    size_w = prod(nw_theta+1) * prod(nw_omega+1);

    variable V(size_v,size_v) hermitian
    variable W(size_w,size_w) hermitian
    variable S0(prod(nS0_theta+ones) * prod(nS0_omega+ones), prod(nS0_theta+ones) * prod(nS0_omega+ones)) hermitian %semidefinite
    variable S1(prod(nS1_theta+ones) * prod(nS1_omega+ones), prod(nS1_theta+ones) * prod(nS1_omega+ones)) hermitian %semidefinite 


    minimize(0)

    subject to
        V >= 0;
        W >= 0;
        S0 >= 0;
        S1 >= 0;

        %V - (0.01) * normPowerGramSparse(decaypower_p,nv_omega,nv_theta) == S0 + selection_matrix(nS1_omega,nS1_theta,nq1_omega,nq1_theta)' * kron(S1,Q1) * selection_matrix(nS1_omega,nS1_theta,nq1_omega,nq1_theta);

        % --- build RHS with V INSIDE (exact) ---
        RHS  = RHSbracketmatrix( ...
                    nv_omega, nv_theta, nf_omega, nf_theta, V, F);

        XRHSX = X' * RHS * X;

        XS1Q1=selection_matrix(nS1_omega,nS1_theta,nq1_omega,nq1_theta);
        
        Gamma=normPowerGramSparse(decaypower_p,nv_omega,nv_theta);

        % --- trace constraints ---
        nConstrv
        nConstr
        for idx = 1:nConstr
            trace( ConstrData{idx}.H * (XRHSX - W) ) == 0;
        end

        for idx = 1:nConstrv
            trace( ConstrDatav{idx}.H * (V - (0.01) * Gamma - S0 - XS1Q1' * kron(S1,Q1) * XS1Q1) )==0;
        end

        trace(V) == 23;
        %V(prod(nv_theta+ones) * (prod(nv_omega+ones)-1) +1:prod(nv_theta+ones) * prod(nv_omega+ones) , prod(nv_theta+ones) * (prod(nv_omega+ones)-1) + 1 :prod(nv_theta+ones) * prod(nv_omega+ones))>=0.0001 * eye(prod(nv_theta+ones));
        %trace(V(prod(nv_theta+ones) * (prod(nv_omega+ones)-1) +1:prod(nv_theta+ones) * prod(nv_omega+ones) , prod(nv_theta+ones) * (prod(nv_omega+ones)-1) + 1 :prod(nv_theta+ones) * prod(nv_omega+ones)))>=0.0001;
 cvx_end

% ================= OUTPUT =================
status_main = cvx_status;
Vmain = V;
Wmain = W;

end

% ===== LOCAL FUNCTIONS (VISIBLE TO WORKERS) =====

function data = curlyH_operator(nr_omega, nr_theta, eta, k)
    c = length(nr_omega);
    d = length(nr_theta);
    H = speye(1);

    for j = 1:d
        H = kron(toeplitz_matrix(nr_theta(j)+1, k(j)), H);
    end

    for i = 1:c
        H = kron(hankel_matrix(nr_omega(i)+1, eta(i)), H);
    end

    data.H = H;
end

function T = toeplitz_matrix(n, k)
    if abs(k) > n-1
        T = sparse(n,n);
        return
    end
    T = spdiags(ones(n,1), k, n, n);
end

function H = hankel_matrix(n, eta)
    if eta < 0 || eta > 2*(n-1)
        H = sparse(n,n);
        return
    end
    k = eta - (n-1);
    H = flipud(spdiags(ones(n,1), k, n, n));
end
