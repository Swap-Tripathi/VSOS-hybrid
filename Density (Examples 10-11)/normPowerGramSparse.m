function P = normPowerGramSparse(p, dimsX1Omega, dimsX1Theta)
% normPowerGramSparse
%
% Sparse diagonal Gram matrix for
%
%   ||omega||^(2p) = (omega_1^2 + ... + omega_c^2)^p
%
% with respect to the standard hybrid basis in equation (3).
%
% INPUT
%   p          : nonnegative integer scalar
%   nv_omega   : vector of maximum omega powers in psi
%   nv_theta   : vector of maximum theta Fourier powers in psi
%
% OUTPUT
%   P          : sparse PSD Gram matrix
%
% The construction uses the multinomial expansion:
%
%   (sum_j omega_j^2)^p
%     = sum_{|alpha|=p} multinomial(p,alpha)*(omega^alpha)^2
%
% Hence P is diagonal.

    % ---------- Input checks ----------
    if ~isscalar(p) || p < 0 || p ~= floor(p)
        error('p must be a nonnegative integer scalar.');
    end

    dimsX1Omega = dimsX1Omega(:).';
    dimsX1Theta = dimsX1Theta(:).';

    c = length(dimsX1Omega);

    % ---------- Multi-indices |alpha| = p ----------
    if c == 1
        alpha = p;
    else
        % Stars-and-bars construction
        bars = nchoosek(1:(p+c-1), c-1);

        alpha = diff([zeros(size(bars,1),1), ...
                      bars, ...
                      (p+c)*ones(size(bars,1),1)], 1, 2) - 1;
    end

    % Check that psi contains every required omega^alpha
    if any(max(alpha,[],1) > dimsX1Omega)
        error(['nv_omega is too small: psi must contain all ', ...
               'omega^alpha with sum(alpha) = p.']);
    end

    % ---------- Size of hybrid basis ----------
    dims = [dimsX1Theta + 1, dimsX1Omega + 1];

    N = prod(dims);

    % MATLAB/Kronecker ordering:
    % first listed dimension varies fastest
    stride = [1, cumprod(dims(1:end-1))];

    nTerms = size(alpha,1);

    idx = zeros(nTerms,1);
    val = zeros(nTerms,1);

    % ---------- Build nonzero diagonal entries ----------
    for r = 1:nTerms

        % All theta powers are zero.
        % Omega powers are alpha(r,:).
        subs = zeros(1,length(dims));
        subs(length(dimsX1Theta) + (1:c)) = alpha(r,:);

        % Linear index of omega^alpha in psi
        idx(r) = 1 + sum(subs .* stride);

        % Multinomial coefficient:
        % p!/(alpha_1! ... alpha_c!)
        val(r) = factorial(p) / prod(factorial(alpha(r,:)));

    end

    % ---------- Sparse diagonal Gram matrix ----------
    P = sparse(idx, idx, val, N, N);

end