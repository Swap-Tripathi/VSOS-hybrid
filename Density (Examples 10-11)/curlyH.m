function data = curlyH(nr_omega, nr_theta, eta, k)
% curlyH
% -------------------------------------------------------------------------
% Builds the numeric sparse operator H(eta,k) such that
%
%   trace( H(eta,k) * R )
%
% equals the original curlyH(R,...) definition.
%
% This function is:
%   - NUMERIC ONLY
%   - PARFOR SAFE
%   - CVX FREE
%
% Output:
%   data.H : sparse matrix
% -------------------------------------------------------------------------

    % Dimensions
    c = length(nr_omega);    % omega dimension
    d = length(nr_theta);    % theta dimension

    % Start with scalar identity
    H = speye(1);

    % ---------------- Toeplitz blocks (theta) ----------------
    for j = 1:d
        Tj = toeplitz_matrix(nr_theta(j)+1, k(j));
        H  = kron(Tj, H);
    end

    % ---------------- Hankel blocks (omega) ----------------
    if c~=0
        for i = 1:c
            Hi = hankel_matrix(nr_omega(i)+1, eta(i));
            H  = kron(Hi, H);
        end
    end
    % Store operator
    data.H = H;
end

% ========================================================================
% Helper: Toeplitz matrix
% ========================================================================
function T = toeplitz_matrix(n, k)

    if abs(k) > n-1
        T = sparse(n,n);
        return
    end

    v = ones(n,1);
    T = spdiags(v, k, n, n);
end

% ========================================================================
% Helper: Hankel matrix
% ========================================================================
function H = hankel_matrix(n, eta)

    if eta < 0 || eta > 2*(n-1)
        H = sparse(n,n);
        return
    end

    k = eta - (n-1);
    v = ones(n,1);
    T = spdiags(v, k, n, n);
    H = flipud(T);
end
