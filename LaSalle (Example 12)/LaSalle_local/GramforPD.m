function V = GramforPD(nv_omega, nv_theta)
% standardVGram
%
% Returns a sparse Hermitian Gram matrix V for
%
%   V(omega,theta) = ||omega||^2 ...
%                  + sum_j (1 - cos(theta_j))
%
% with respect to the standard hybrid basis in equation (3).
%
% INPUTS
%   nv_omega : 1 x c vector of Euclidean representation degrees
%   nv_theta : 1 x d vector of trigonometric representation degrees
%
% OUTPUT
%   V : sparse Hermitian Gram matrix satisfying
%
%       psi' * V * psi
%         = sum_j omega_j^2
%           + sum_j (1 - cos(theta_j))
%
% REQUIREMENTS
%   nv_omega(j) >= 1 for every omega variable
%   nv_theta(j) >= 1 for every theta variable

    nv_omega = nv_omega(:).';
    nv_theta = nv_theta(:).';

    c = length(nv_omega);
    d = length(nv_theta);

    % Check that required monomials exist in the basis
    if any(nv_omega < 1)
        error('All entries of nv_omega must be at least 1.');
    end

    if any(nv_theta < 1)
        error('All entries of nv_theta must be at least 1.');
    end

    % ---------------------------------------------------------
    % Basis dimensions
    %
    % Ordering consistent with equation (3):
    %
    % theta_1 varies fastest, followed by theta_2, ...
    % then omega_1, omega_2, ...
    % ---------------------------------------------------------

    dims = [nv_theta + 1, nv_omega + 1];

    N = prod(dims);

    stride = [1, cumprod(dims(1:end-1))];

    % ---------------------------------------------------------
    % Preallocate sparse matrix
    %
    % Nonzeros:
    %   1 constant entry
    %   c diagonal omega entries
    %   2d off-diagonal theta entries
    % ---------------------------------------------------------

    V = sparse(N,N);

    % ---------------------------------------------------------
    % Constant term:
    %
    % sum_j 1 = d
    % ---------------------------------------------------------

    V(1,1) = d;

    % ---------------------------------------------------------
    % ||omega||^2
    %
    % Each omega_j^2 is obtained from the diagonal Gram
    % entry corresponding to the basis monomial omega_j.
    % ---------------------------------------------------------

    for j = 1:c

        subs = zeros(1,c+d);

        % exponent 1 on omega_j
        subs(d+j) = 1;

        idx = 1 + sum(subs .* stride);

        V(idx,idx) = V(idx,idx) + 1;

    end

    % ---------------------------------------------------------
    % -sum_j cos(theta_j)
    %
    % Since
    %
    %   cos(theta_j)
    %      = 0.5*exp(i*theta_j)
    %      + 0.5*exp(-i*theta_j),
    %
    % place -1/2 between the constant basis element and
    % exp(i*theta_j).
    % ---------------------------------------------------------

    for j = 1:d

        subs = zeros(1,c+d);

        % Fourier exponent 1 for theta_j
        subs(j) = 1;

        idx = 1 + sum(subs .* stride);

        V(1,idx) = V(1,idx) - 0.5;
        V(idx,1) = V(idx,1) - 0.5;

    end

end