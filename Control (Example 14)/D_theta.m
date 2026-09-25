function D = D_theta(l, n_omega, n_theta)

c = numel(n_omega);
d = numel(n_theta);

D = speye(1);


for j = 1:d
    if j == l
        D = kron(1i * diag(0:n_theta(j)), D);
    else
        D = kron(speye(n_theta(j)+1), D);
    end
end

% omega block
for j = 1:c
    D = kron(speye(n_omega(j)+1), D);
end

end