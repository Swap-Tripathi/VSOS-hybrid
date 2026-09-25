%syms omega c s; createMonomialBasis(1,3, omega, c, s)
% SOSDEMO --- Lyapunov Density Search after lifting


clear; clc;
tic;
sostools_root = 'ADD ROOT HERE';
addpath(genpath(sostools_root));
verbose = 1;
nv_omega=5;
nv_theta=3;
epsilon=0.01;
R=4;
pvar omega theta c s;
vars = [omega; theta; c; s];

% Constructing the vector field dx/dt = f
f = [-0.8*omega-10*s;
     omega;
    -s*omega;
    c*omega];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);

% =============================================
% The Lyapunov function V(x):
Zv=createMonomialBasis(nv_omega, nv_theta, omega, c, s);
[prog,V] = sospolyvar(prog,Zv);

% The polynomial lambdaV:
Zlv=createMonomialBasis(nv_omega, max(nv_theta-2,0), omega, c, s);
[prog,lv] = sospolyvar(prog,Zlv);

nw_omega=nv_omega+1;
nw_theta=nv_theta+1;

% The Lyapunov function W(x):
Zw=createMonomialBasis(nw_omega, nw_theta, omega, c, s);
[prog,W] = sospolyvar(prog,Zw);

% The polynomial lambdaW:
Zlw=createMonomialBasis(nw_omega, max(nw_theta-2,0), omega, c, s);
[prog,lw] = sospolyvar(prog,Zlw);

% The polynomial lambdaVtilde:
Zlvtilde = createMonomialBasis(nv_omega, max(2*ceil(nv_theta/2)-2,0), omega,c,s);
[prog,lvtilde] = sospolyvar(prog,Zlvtilde);

% The polynomial S1V:
%Zs1v=createMonomialBasis(nv_omega-1, nv_theta, omega, c, s);
%[prog,s1v] = sossosvar(prog,Zs1v);
Zs1v = createSOSBasis(nv_omega-1, ceil(nv_theta/2), omega, c, s);
[prog,s1v] = sossosvar(prog,Zs1v);

% The polynomial lambdaDensity:
[prog,ld] = sospolyvar(prog,Zlw);

% =============================================
% Next, define SOSP constraints

% Constraint 1 : V(x) - lambdaV*(c^2+s^2-1) >= 0
prog = sosineq(prog,V-lv*(c^2+s^2-1));

% Constraint 2 : W(x) - lambdaW*(c^2+s^2-1) >= 0
prog = sosineq(prog,W-lw*(c^2+s^2-1));

% Constraint 3 : V(x) >= s1v (omega^2-R^2)+epsilon omega^2+ lvtilde*(c^2+s^2-1)
prog = sosineq(prog,V-s1v*(omega^2-R^2)-epsilon*omega^2-lvtilde*(c^2+s^2-1));

% Constraint 4: W = -grad(V)*f + V*div(f)
expr = -(diff(V,omega)*f(1)+diff(V,theta)*f(2)+diff(V,c)*f(3)+diff(V,s)*f(4))+V*(diff(f(1),omega)+diff(f(3),c)+diff(f(4),s))-W- ld*(c^2+s^2-1);
prog = soseq(prog,expr);

Veq1 = subs(V,omega,0);
Veq1 = subs(Veq1,c,1);
Veq1 = subs(Veq1,s,0);

Veq2 = subs(V,omega,0);
Veq2 = subs(Veq2,c,-1);
Veq2 = subs(Veq2,s,0);

prog = soseq(prog,Veq1);
prog = soseq(prog,Veq2);

% =============================================
% And call solver
solver_opt.solver = 'mosek';
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-6;    % precision low equivalent 
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-6;    % precision low equivalent
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6;    % precision low equivalent
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-6;    % precision low equivalent

prog = sossolve(prog,solver_opt);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,V);

total_time = toc;   % Stop timer

% ============================================================
% Recover V and convert c=cos(theta), s=sin(theta)
% ============================================================

% Convert V to symbolic form and substitute c=cos(theta), s=sin(theta)
% (sym() cannot directly convert a SOSTOOLS 'polynomial' object, so we
% rebuild it term-by-term from its coefficient/degmat data instead.)
%
% IMPORTANT: SOSTOOLS/multipoly 'polynomial' objects store variables in
% ALPHABETICAL order internally (regardless of the order you declared
% them with `pvar omega c s`). Alphabetically: c < omega < s. So
% SOLV.degmat columns are ordered [c, omega, s], NOT [omega, c, s].
% We must look up the correct column for each variable via SOLV.varname
% instead of assuming a fixed column order -- using the wrong order is
% what caused the omega-exponent to leak into the cos/sin powers
% (e.g. cos(10*theta) artifacts).

syms Om theta real
tol = 1e-8;

idx_omega = find(strcmp(SOLV.varname,'omega'));
idx_c     = find(strcmp(SOLV.varname,'c'));
idx_s     = find(strcmp(SOLV.varname,'s'));

Vsym = sym(0);
for ii = 1:length(SOLV.coefficient)
    coeff  = full(SOLV.coefficient(ii));
    powers = full(SOLV.degmat(ii,:));
    if abs(coeff) > tol
        po = powers(idx_omega);
        pc = powers(idx_c);
        ps = powers(idx_s);
        Vsym = Vsym + coeff * Om^po * cos(theta)^pc * sin(theta)^ps;
    end
end

Vtheta = simplify(Vsym, 'Steps', 100);

fprintf('\n=============================================\n');
fprintf('V(omega,theta) after c=cos(theta), s=sin(theta):\n');
fprintf('=============================================\n');
%disp(Vtheta)
vpa(Vtheta,6)

% Plot 1/V(omega,theta) as an "energy" surface
Efun = matlabFunction(log(abs(1/Vtheta)), 'Vars', [Om, theta]);

om_range    = linspace(-4, 4, 200);
theta_range = linspace(-pi, pi, 200);
[OM, TH]    = meshgrid(om_range, theta_range);

E = Efun(OM, TH);

figure;
pcolor(OM, TH, E);
shading interp;
xlabel('\omega');
ylabel('\theta');
title('Energy plot: 1/V(\omega,\theta)');
colorbar;

function Z = createMonomialBasis(n_omega,n_theta,omega,c,s)

    Z = [];

    for i = 0:2*n_omega
        for j = 0:n_theta
            for k = 0:(n_theta-j)
                Z = [Z; omega^i*c^j*s^k];
            end
        end
    end
end

function Z = createSOSBasis(n_omega,n_theta,omega,c,s)
% Creates a Gram monomial basis for an SOS polynomial.
%
% Resulting SOS degree:
%   deg_omega <= 2*n_omega
%   deg_(c,s)  <= 2*n_theta
%
% n_omega and n_theta are the degrees of the Gram vector.

    Z = [];

    for i = 0:n_omega
        for j = 0:n_theta
            for k = 0:(n_theta-j)
                Z = [Z; omega^i*c^j*s^k];
            end
        end
    end
end
