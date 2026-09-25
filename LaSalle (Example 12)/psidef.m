%uncomment to check the validity of the function
%syms omega1 omega2 theta1 theta2 real;
%psidef([omega1,omega2],[theta1,theta2],[1 2],[3,4])

function psidef = psidef(omega, theta, n_omega, n_theta)
    
    c = numel(n_omega);
    d = numel(n_theta);
    % z[Theta] = [exp(i*k*theta), k = 0..degTrig]
    psidef=1; 

    for i=1:d
    psidef = kron( exp(1i * (0:n_theta(i)) * theta(i)), psidef);
    end
    % y[omega] = [omega^k, k = 0..degPoly]
    for i=1:c
    psidef = kron( omega(i) .^ (0:n_omega(i)) , psidef);
    end
    % Flatten into a row vector
    psidef = psidef(:);   % transpose (') to make row vector

end