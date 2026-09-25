clc; clear; close; close all; clear functions;
cvx_clear;

% example two second-order coupled oscillators

example=14;                          % change example number

nf_omega = [1];                     % change according to system (representation-size vector)
nf_theta = [1];                     % change according to system (representation-size vector)

c = numel(nf_omega);
d = numel(nf_theta);

F = cell(1,c+d);

% Call F{1},...F{c+d} for the examples in the paper

if example==14
        F{1}= sparse([1     1       2       3   ], ...        % row indices
                     [2     3       1       1   ], ...        % column indices
                     [5*1j  -0.4    -5*1j   -0.4], ... % values
                     prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        F{2}= sparse([3 1], [1 3], [1/2 1/2], prod(nf_theta+ones) * prod(nf_omega+ones), prod(nf_theta+ones) * prod(nf_omega+ones));
        nv_omega_common=1;
        nv_theta_common=1;
end

nv_omega = [1] * nv_omega_common;     % Start with the smallest value and increase if infeasible
nv_theta = [1] * nv_theta_common;     % Start with the smallest value and increase if infeasible


    disp(['============ START (',datestr(now),') ============']);
    tic

    [status_main, Vgram, Wgram] = SDPLasalle(nv_omega, nv_theta, nf_omega, nf_theta, F);
    disp(status_main);
    toc

    syms omega theta real; 
    expr=psidef(omega,theta,1,1)' * Vgram * psidef(omega,theta,1,1);
    expr = expand(expr);
    expr = simplify(vpa(rewrite(expr, 'sincos'),6))