function [status_main, Vmain, Wmain, S0W, S1W, S2W] = SDPLasalle_local( ...
        nv_omega, nv_theta, nf_omega, nf_theta, nq_omega, nq_theta, F,Q)

% ================= DIMENSIONS =================
nw_omega = nv_omega + nf_omega;
nw_theta = nv_theta + nf_theta;

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

% ================= SELECTION MATRIX =================
X = selection_matrix(nv_omega, nv_theta, nf_omega, nf_theta);

% ================= PRECOMPUTE curlH OPERATORS =================
Klist   = STheta(:,1:(end+1)/2).';
Etalist = SOmega.';

nK   = size(Klist,1);
nEta = size(Etalist,1);
nConstr = nK * nEta;

ConstrData = cell(nConstr,1);

parfor idx = 1:nConstr
    [iK,iEta] = ind2sub([nK,nEta],idx);
    k   = Klist(iK,:).';
    eta = Etalist(iEta,:).';
    ConstrData{idx} = curlyH_operator(nw_omega, nw_theta, eta, k);
end

% ================= SDP =================
cvx_clear
cvx_solver mosek
cvx_precision low

cvx_begin sdp quiet

    size_v = prod(nv_theta+1) * prod(nv_omega+1);
    size_s0w = prod(nw_theta+1) * prod(nw_omega+1);
    size_s12w = prod(nw_theta-nq_theta+1) * prod(nw_omega-nq_omega+1);

    variable V(size_v,size_v) hermitian
    variable S0W(size_s0w,size_s0w) hermitian
    variable S1W(size_s12w,size_s12w) hermitian
    variable S2W(size_s12w,size_s12w) hermitian

    
    minimize(0)

    subject to
        V >= 0;
        S0W >= 0;
        S1W >= 0;
        S2W >= 0;

        % --- build RHS with V INSIDE (exact) ---
        RHS  = RHSbracketmatrix_LaSalle_local( ...
                    nv_omega, nv_theta, nf_omega, nf_theta, V, F, Q);

        XRHSX = X' * RHS * X;
        
        W=S0W+X' * ( kron(S1W , Q{1}) + kron(S2W , Q{2}) ) * X;

        % --- trace constraints ---
        for idx = 1:nConstr
            trace( ConstrData{idx}.H * (XRHSX - (0.01) * GramforPD(nv_omega+nf_omega,nv_theta+nf_theta) - W) ) == 0;
        end

        trace(V) >= 5;
        psidef(0,0,nv_omega+nf_omega,nv_theta+nf_theta)' * S0W * psidef(0,0,nv_omega+nf_omega,nv_theta+nf_theta)==0;
        psidef(0,0,nv_omega+nf_omega-nq_omega,nv_theta+nf_theta-nq_theta)' * S1W * psidef(0,0,nv_omega+nf_omega-nq_omega,nv_theta+nf_theta-nq_theta)==0;
        psidef(0,0,nv_omega+nf_omega-nq_omega,nv_theta+nf_theta-nq_theta)' * S2W * psidef(0,0,nv_omega+nf_omega-nq_omega,nv_theta+nf_theta-nq_theta)==0;
cvx_end

% ================= OUTPUT =================
status_main = cvx_status;
Vmain = V;
Wmain = S0W+X' * ( kron(S1W , Q{1}) + kron(S2W , Q{2}) ) * X;
S0W = S0W;
S1W = S1W;
S2W = S2W;
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
