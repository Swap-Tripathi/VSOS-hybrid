function [status_main, Vmain, Wmain] = SDPmain(nutilde_omega, nutilde_theta, ns_omega, ns_theta, S, R, example)
    
    nw_omega =  nutilde_omega + ns_omega;
    nw_theta =  nutilde_theta + ns_theta;
    c = numel(nw_omega);
    d = numel(nw_theta);
    
    %--------GENERATE ALL INDICES BETWEEN -nw_theta and nw_theta-------------
    grids = cell(1,d);
        for j = 1:d
            grids{j} = -nw_theta(j) : nw_theta(j);
        end
    [Kt{1:d}] = ndgrid(grids{:});
    STheta = zeros(numel(Kt{1}), d);
        for j = 1:d
            STheta(:,j) = Kt{j}(:);
        end
    STheta=STheta.';

    %--------GENERATE ALL INDICES BETWEEN 0 and 2nw_omega-------------
    grids = cell(1,c);
        for j = 1:c
            grids{j} = 0 : 2*nw_omega(j);
        end
    [Ko{1:c}] = ndgrid(grids{:});
    SOmega = zeros(numel(Ko{1}), c);
        for j = 1:c
            SOmega(:,j) = Ko{j}(:);
        end
    SOmega=SOmega.';
    
    X = selection_matrix(ns_omega,ns_theta,nutilde_omega,nutilde_theta); % ordering of s and utilde in the expression should match with ordering here
                                                                         % if S comes to the left of Utilde in
                                                                         % tensor product in the theorem,
                                                                         % then n_s comes before n_utilde in X.
    
    if example==16                                                        % bounds will be used in the SDP to make SDP matrices sufficiently far away from 0 matrix.
        bound=10;
    elseif example==17
        bound=2595;
    end

    % -----------SDP starts here------------
    cvx_clear
    cvx_solver mosek % can change to sedumi or sdpt3: results vary
    cvx_precision high
    cvx_begin sdp quiet
    
    size_utilde = prod(nutilde_theta+ones) * prod(nutilde_omega+ones);
    size_w = prod(nw_theta+ones) * prod(nw_omega+ones);
    variable Utilde(size_utilde, size_utilde) hermitian %semidefinite
    variable Wtilde(size_w, size_w) hermitian %semidefinite  
    
    minimize(0);
    subject to
        (Utilde);
        (Wtilde)>=0;
        
        for k = STheta(:,1:(end+1)/2)

            for eta = SOmega(:,:)
                
                curlyH( X' * RHSbracketmatrix(nutilde_omega,nutilde_theta,ns_omega,ns_theta, S, R, Utilde)  * X - Wtilde , nw_omega, nw_theta , eta , k)==0;
            
            end

        end
        trace(Wtilde) >= bound;
        %psidef(0,0, nv_omega, nv_theta)' * V * psidef(0,0, nv_omega, nv_theta) ==0;
    cvx_end

    status_main = cvx_status;
    Vmain = Utilde;%round(full(V),3);    % converts sparse matrix to a normal one
    Wmain = Wtilde;%round(full(W),3);    % converts sparse matrix to a normal one
end
