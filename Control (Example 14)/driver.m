clc; clear; close; close all; clear functions;
cvx_clear;

% Swing-up control for pendulum on cart

example=14;                          % 16 for pendulum with fixed cart, 17 for moving cart.

ns_omega = [1];                     % change according to system (representation-size vector)
ns_theta = [2];                     % change according to system (representation-size vector)
% n_r=n_r imposed later

% define systems

if example==16                       % pendulum with fixed cart
        a=1;
        S = sparse([1       1     1       2       3     4       5], ...        % row indices
                   [1       2     3       1       1     5       4], ...        % column indices
                   [1/2    -1/2   1/4    -1/2     1/4   1/4     1/4], ... % values
                     prod(ns_theta+ones) * prod(ns_omega+ones), prod(ns_theta+ones) * prod(ns_omega+ones));
        R = sparse([1 5], ...
                   [5 1], ...
                   [-1/2 -1/2], prod(ns_theta+ones) * prod(ns_omega+ones), prod(ns_theta+ones) * prod(ns_omega+ones)); % n_r=n_s
        nutilde_omega_common=1;
        nutilde_theta_common=1;
elseif example==14
        a=185.094;
        S = sparse([1       1     1       2       3     4       5], ...        % row indices
                   [1       2     3       1       1     5       4], ...        % column indices
                   [-a/2    a/2   -a/4    a/2     -a/4  -1/4    -1/4], ... % values
                     prod(ns_theta+ones) * prod(ns_omega+ones), prod(ns_theta+ones) * prod(ns_omega+ones));
        R = sparse([1   5], ...
                   [5   1], ...
                   [1/2 1/2], prod(ns_theta+ones) * prod(ns_omega+ones), prod(ns_theta+ones) * prod(ns_omega+ones)); % n_r=n_s
        nutilde_omega_common=1;
        nutilde_theta_common=1;
end

nutilde_omega = [1] * nutilde_omega_common;     % Start with the smallest value and increase if infeasible
nutilde_theta = [1] * nutilde_theta_common;     % Start with the smallest value and increase if infeasible


    disp(['============ START (',datestr(now),') ============']);
    tic

    [status_main, Utildegram, Wtildegram] = SDPmain(nutilde_omega, nutilde_theta, ns_omega, ns_theta, S, R, example); % n_s=n_r
    disp(status_main);
    toc
    
    %% ========== Run separately later for PLOTTING log(1/V) =============
    
    if example==16
        omega_axis=7;
        density_streamlines=50;
    elseif example==17
        omega_axis=60;
        density_streamlines=1000;
    end
    syms omega1 theta1 real;
    Utilde_expr=simplify(rewrite(psidef([omega1],[theta1],nutilde_omega,nutilde_theta)' * round(Utildegram,5) * psidef([omega1],[theta1],nutilde_omega,nutilde_theta),'sincos'));
    if example==16
        k_expr=cos(theta1);
    elseif example==17
        k_expr=-cos(theta1);
    end
    F1expr = a*sin(theta1)+(0.5*omega1^2+a*(cos(theta1)-1))*Utilde_expr*k_expr;
    F2expr = omega1;
    Vexpr  = (0.5*omega1^2+a*(cos(theta1)-1))^2;

    if strcmp(status_main,'Solved')
        V_num = matlabFunction(Vexpr, 'Vars', [theta1 omega1]);
        F1num = matlabFunction(F1expr, 'Vars', [theta1 omega1]);
        F2num = matlabFunction(F2expr, 'Vars', [theta1 omega1]);
        theta_vals = linspace(-pi, pi, 400);
        omega_vals = linspace(-omega_axis, omega_axis, 400);
        
        [TH, OM] = meshgrid(theta_vals, omega_vals);
        
        Z = abs(V_num(TH, OM));
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
        
        xlabel('\theta')
        ylabel('\omega')
        title('Energy map of log \rho(\omega,\theta) with vector field')
        
        axis tight
     
    end

    %disp(['============= END (',datestr(now),') ============='])
















