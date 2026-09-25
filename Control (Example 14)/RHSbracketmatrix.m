function RHSbracketmatrix = RHSbracketmatrix(nutilde_omega, nutilde_theta, ns_omega, ns_theta, S, R, Utilde)

rsize = prod(nutilde_theta+ones) * prod(nutilde_omega+ones) * prod(ns_theta+ones) * prod(ns_omega+ones);

RHSbracketmatrix=sparse(rsize, rsize);
    
Dw_utilde = D_omega(1,nutilde_omega,nutilde_theta);
RHSbracketmatrix = RHSbracketmatrix + kron(S , Dw_utilde' * Utilde + Utilde * Dw_utilde)+ kron(R,Utilde);
end