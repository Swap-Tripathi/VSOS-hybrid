clc; clear; close; close all; clear functions;
cvx_clear;


% example two second-order coupled oscillators

example=10;                         % change example number (11, 12 or 13)

R=4;

nf_omega = [1];                     % change according to system (representation-size vector)
nf_theta = [1];                     % change according to system (representation-size vector)

c = numel(nf_omega);
d = numel(nf_theta);

decaypower_p=ceil((length(nf_omega)+1)/2); % use while plugging in normPowerGramSparse. this is the smallest p implying integrability of v(omega,theta)>=epsilon*||omega||^{2p}

F = cell(1,c+d);

% Call F{1},...F{c+d} for the examples in the paper

if example==10
        F{1}= sparse([1     1       2       3   ], ...        % row indices
                     [2     3       1       1   ], ...        % column indices
                     [5*1j  -0.4    -5*1j   -0.4], ... % values
                     prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        F{2}= sparse([3 1], [1 3], [1/2 1/2], prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        nv_omega_common=5; % make sure polypower>=decaypower_p;
        nv_theta_common=3;
        nv_omega = [1] * nv_omega_common;     % Start with the smallest value and increase if infeasible
        nv_theta = [1] * nv_theta_common;     % Start with the smallest value and increase if infeasible
elseif example==11
        Kosc= 10;   
        alpha_osc=5; 
        om12= 10;    
        F{1}= sparse([1         1           1                   2           3   ], ...        % row indices
                     [1         2           3                   1           1   ], ...        % column indices
                     [om12      Kosc*1j/2   -alpha_osc/2      -Kosc*1j/2    -alpha_osc/2], ... % values
                     prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        F{2}= sparse([3 1], [1 3], [1/2 1/2], prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        nv_omega_common=5; % make sure polypower>=decaypower_p;
        nv_theta_common=3;
        nv_omega = [1] * nv_omega_common;     % Start with the smallest value and increase if infeasible
        nv_theta = [1] * nv_theta_common;     % Start with the smallest value and increase if infeasible
end

for i=1:1:c+d
    if isequal(size(F{i}), [0 0])
    error('==== Error: empty component in the vector field ====')
    return
    end
end



    disp(['============ START (',datestr(now),') ============']);
    tic

    [status_main, Vgram, Wgram] = SDPmain(nv_omega, nv_theta, nf_omega, nf_theta, F, decaypower_p,R);
    disp(status_main);
    toc
    
    %% ========== Run separately later for PLOTTING log(1/V) =============
    
    omega_axis=10;
    density_streamlines=10;
    
    syms omega1 theta1 real; % change variables dependent on the example number. Remove omegas for pure rotational systems (example 11)
    Vexpr=simplify(rewrite(psidef([omega1],[theta1],nv_omega,nv_theta)' * Vgram * psidef([omega1],[theta1],nv_omega,nv_theta),'sincos'));
    F1expr = rewrite(psidef([omega1],[theta1],nf_omega,nf_theta)' * F{1} * psidef([omega1],[theta1],nf_omega,nf_theta),'sincos');
    F2expr = rewrite(psidef([omega1],[theta1],nf_omega,nf_theta)' * F{2} * psidef([omega1],[theta1],nf_omega,nf_theta),'sincos');


    if strcmp(status_main,'Solved')
        Vnum = matlabFunction(Vexpr, 'Vars', [theta1 omega1]);
        F1num = matlabFunction(F1expr, 'Vars', [theta1 omega1]);
        F2num = matlabFunction(F2expr, 'Vars', [theta1 omega1]);
        theta_vals = linspace(-pi, pi, 400);
        omega_vals = linspace(-omega_axis, omega_axis, 400);
        
        [TH, OM] = meshgrid(theta_vals, omega_vals);
        
        Z = abs(Vnum(TH, OM));
        Z(Z < 1e-12) = 1e-12;     % regularization
        E = -log(Z);              % log(1/V)
        % VECTOR FIELD GRID (coarse)
        theta_q = linspace(-pi, pi, density_streamlines);
        omega_q = linspace(-omega_axis, omega_axis, density_streamlines);
        [THq, OMq] = meshgrid(theta_q, omega_q);
        
        U = F2num(THq, OMq);
        V = F1num(THq, OMq);
        
        % PLOT
        figure
        imagesc(theta_vals, omega_vals, E)
        set(gca,'YDir','normal')
        hold on
        colorbar
        
        streamslice(THq, OMq, U, V);
        set(streamslice(THq, OMq, U, V), 'Color', 'k', 'LineWidth', 0.8)
        %quiver(THq, OMq, U, V, 'k', 'LineWidth', 1)
        
        xlabel('\theta')
        ylabel('\omega')
        title('Energy map of log \rho(\omega,\theta) with vector field')
        
        axis tight
     
    end

    %disp(['============= END (',datestr(now),') ============='])
















