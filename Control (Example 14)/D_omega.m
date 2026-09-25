function D = D_omega(l, n_omega, n_theta)

c = numel(n_omega);
d = numel(n_theta);

D = speye(1);


% theta block (all identities)
for j = 1:d
    D = kron(speye(n_theta(j)+1), D);
end

% omega block
for j = 1:c
    if j == l
        D = kron(diag(1:n_omega(j), -1), D);
    else
        D = kron(speye(n_omega(j)+1), D);
    end
end

end