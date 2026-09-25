function X = selection_matrix(nvomega,nvtheta, npomega,nptheta)

Ev  = exponent_table(nvomega , nvtheta);      % order changed due to how ndgrid works. The first grid (n_{v_theta}(1)) due to the change changes the quickest
Ep  = exponent_table(npomega , nptheta);      % |psi_nv| x (c+d) where each row is a multiindex ordered by tensor product, exponent of theta_1 changes quickest
Evp = exponent_table(nvomega + npomega , nvtheta + nptheta);

Nv  = size(Ev,1);                              % fetches the size of the psi_{n_v}  
Np  = size(Ep,1);                              % fetches the size of the psi_{n_p}  
Nvp = size(Evp,1);                             % fetches the size of the psi_{n_v+n_p}  

X = sparse(Nv*Np, Nvp);                        % sparse only saves positions of nonzero values and saves computation time.

row = 1;
for i = 1:Nv
    for j = 1:Np
        exp_sum = Ev(i,:) + Ep(j,:);
        k = find(all(Evp == exp_sum, 2), 1);
        X(row, k) = 1;
        row = row + 1;
    end
end

end

function E = exponent_table(nomega,ntheta) % produces the standard basis vectors as grids
% n : 1 x (c+d) vector of nonnegative integers
% E : prod(vec+1) x (c+d) exponent table (Kronecker order). These are all possible exponents...
% (k1,...kd,eta1,...etac) ordered in a way that k1 changes quickest followed by k2 and so on, if...
% the input vector is (n_{v_theta}(1),...,n_{v_theta}(d),n_{v_omega}(1),...,n_{v_omega}(c))

nomega = nomega(:).';                             % ensure row vector                
ntheta = ntheta(:).';                             % ensure row vector     
c=numel(nomega);                                  % numel(n) is the size of row vector n.
d=numel(ntheta);                                  % numel(n) is the size of row vector n.

grids = cell(1,c+d);               
for k = 1:d
    grids{k} = 0:ntheta(k);
end
for k = 1:c
    grids{d+k} = 0:nomega(k);
end

[G{1:c+d}] = ndgrid(grids{:});       % grids{:} lists all cell array as comma separated
                                        % ndgrid(grids{1}, grids{2}, ..., grids{c+d}) construct
                                        % (c+d)-dimensional coordinate grids with all values an
                                        % exponent corresponding to a variable can take
                                        
E = zeros(numel(G{1}), c+d);         % numel(G{1}) gives total number of combinations. E stores one multiindex per row

for k = 1:c+d
    E(:,k) = G{k}(:);                   % flattens each grid into a column vector.
end
end
