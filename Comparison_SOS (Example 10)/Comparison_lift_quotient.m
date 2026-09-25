% ============================================================
% SOSDEMO --- Lyapunov Density Search after quotient lifting
%
% Variables:
%       omega = angular velocity
%       c     = cos(theta)
%       s     = sin(theta)
%
% Quotient relation:
%       c^2 + s^2 - 1 = 0
%
% Canonical quotient basis:
%       omega^i s^k
%       omega^i c s^k
%
% because
%       c^2 = 1 - s^2.
% ============================================================

clear; clc;
tic;

% ============================================================
% PARAMETERS
% ============================================================

nv_omega = 10;
nv_theta =6;

epsilon = 0.01;
R = 4;

% ============================================================
% SOSTOOLS
% ============================================================

sostools_root = ...
    'ADD ROOT HERE';

addpath(genpath(sostools_root));

verbose = 1;


% ============================================================
% POLYNOMIAL VARIABLES
% ============================================================

pvar omega c s;

vars = [omega; c; s];

% Quotient relation
g = c^2 + s^2 - 1;

% ============================================================
% VECTOR FIELD
%
% dx/dt = f
% ============================================================

f = [ -0.8*omega - 10*s;
      omega;
      -s*omega;
       c*omega ];


% This should return zero.

% ============================================================
% INITIALIZE SOS PROGRAM
% ============================================================

prog = sosprogram(vars);


% ============================================================
% 1. LYAPUNOV FUNCTION V
%
% V belongs to the quotient-space basis
%
% omega^i s^k
% omega^i c s^k
%
% i = 0,...,2*nv_omega
% ============================================================

Zv = createQuotientBasis( ...
    nv_omega, ...
    nv_theta, ...
    omega, c, s);

[prog,V] = sospolyvar(prog,Zv);

% ============================================================
% 2. MULTIPLIER lambdaV
%
% V - lambdaV*g >= 0
% ============================================================

Zlv = createQuotientBasis( ...
    nv_omega, ...
    max(nv_theta-2,0), ...
    omega,c,s);

[prog,lv] = sospolyvar(prog,Zlv);

% ============================================================
% 3. W
%
% W has one additional degree in omega and theta
% ============================================================

nw_omega = nv_omega + 1;
nw_theta = nv_theta + 1;

Zw = createQuotientBasis( ...
    nw_omega, ...
    nw_theta, ...
    omega,c,s);

[prog,W] = sospolyvar(prog,Zw);

% ============================================================
% 4. MULTIPLIER lambdaW
%
% W - lambdaW*g >= 0
% ============================================================

Zlw = createQuotientBasis( ...
    nw_omega, ...
    max(nw_theta-2,0), ...
    omega,c,s);

[prog,lw] = sospolyvar(prog,Zlw);

% ============================================================
% 5. MULTIPLIER lambdaVtilde
%
% Constraint:
%
% V >= s1v*(omega^2-R^2)
%      + epsilon*omega^2
%      + lambdaVtilde*g
%
% ============================================================

Zlvtilde = createQuotientBasis( ...
    nv_omega, ...
    max(2*ceil(nv_theta/2)-2,0), ...
    omega,c,s);

[prog,lvtilde] = sospolyvar(prog,Zlvtilde);

% ============================================================
% 6. SOS POLYNOMIAL s1v
%
% s1v >= 0
%
% We use a quotient-space Gram basis.
% ============================================================

Zs1v = createQuotientSOSBasis( ...
    nv_omega-1, ...
    ceil(nv_theta/2), ...
    omega,c,s);

[prog,s1v] = sossosvar(prog,Zs1v);

% ============================================================
% 7. lambdaDensity
%
% W = -grad(V)*f + V*div(f)
%
% modulo g = 0
%
% Therefore:
%
% -grad(V)*f + V*div(f) - W
% - lambdaDensity*g = 0
%
% ============================================================

[prog,ld] = sospolyvar(prog,Zlw);

% ============================================================
% SOS CONSTRAINTS
% ============================================================

% ------------------------------------------------------------
% Constraint 1:
%
% V >= 0 on c^2+s^2=1
%
% Certificate:
%
% V - lambdaV*g is SOS
% ------------------------------------------------------------

prog = sosineq(prog, V - lv*g);

% ------------------------------------------------------------
% Constraint 2:
%
% W >= 0 on c^2+s^2=1
%
% ------------------------------------------------------------

prog = sosineq(prog, W - lw*g);

% ------------------------------------------------------------
% Constraint 3:
%
% V >= s1v*(omega^2-R^2)
%       + epsilon*omega^2
%
% on c^2+s^2=1
%
% ------------------------------------------------------------

prog = sosineq( ...
    prog, ...
    V ...
    - s1v*(omega^2-R^2) ...
    - epsilon*omega^2 ...
    - lvtilde*g);

% ------------------------------------------------------------
% Constraint 4:
%
% W =
% -grad(V)*f + V*div(f)
%
% modulo g=0
%
% ------------------------------------------------------------

divf = diff(f(1),omega) ...
     + diff(f(2),c) ...
     + diff(f(3),s);

LieDensity = ...
      -(diff(V,omega)*f(1) ...
      + diff(V,c)*f(2) ...
      + diff(V,s)*f(3)) ...
      + V*divf;

expr = LieDensity - W - ld*g;

prog = soseq(prog,expr);

% ============================================================
% NORMALIZATION OF V
%
% V(0,1,0) = 0
% V(0,-1,0) = 0
%
% These are the two equilibrium points corresponding to
% theta = 0 and theta = pi.
% ============================================================

%Veq1 = subs(V,omega,0);
%Veq1 = subs(Veq1,c,1);
%Veq1 = subs(Veq1,s,0);
%
%Veq2 = subs(V,omega,0);
%Veq2 = subs(Veq2,c,-1);
%Veq2 = subs(Veq2,s,0);
%
%prog = soseq(prog,Veq1);
%prog = soseq(prog,Veq2);

% ============================================================
% SOLVE
% ============================================================

solver_opt.solver = 'mosek';

solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-6;    % precision low equivalent 
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-6;    % precision low equivalent
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6;    % precision low equivalent
solver_opt.params.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-6;    % precision low equivalent

prog = sossolve(prog,solver_opt);

% ============================================================
% RECOVER SOLUTION
% ============================================================

SOLV = sosgetsol(prog,V);

% ============================================================
% TIME
% ============================================================

total_time = toc;
% ============================================================
% SOLVER STATUS
% ============================================================

info = prog.solinfo.info;


% ============================================================
% RECOVER V(omega,theta)
% ============================================================

syms Om theta real

tol = 1e-6;

% SOSTOOLS stores variables alphabetically.
idx_omega = find(strcmp(SOLV.varname,'omega'));
idx_c     = find(strcmp(SOLV.varname,'c'));
idx_s     = find(strcmp(SOLV.varname,'s'));

Vsym = sym(0);
flag=0;

for ii = 1:length(SOLV.coefficient)

    coeff  = full(SOLV.coefficient(ii));
    powers = full(SOLV.degmat(ii,:));

    if abs(coeff) > tol
        
        flag=1;
        po = powers(idx_omega);
        pc = powers(idx_c);
        ps = powers(idx_s);

        Vsym = Vsym ...
            + coeff ...
            * Om^po ...
            * cos(theta)^pc ...
            * sin(theta)^ps;
    end
end

Vtheta = simplify(Vsym,'Steps',100);

fprintf('\n=============================================\n');
fprintf('V(omega,theta):\n');
fprintf('=============================================\n');

if flag==1
vpa(Vtheta,6)
else
  fprintf('zero solution');
end
% Plot
% ============================================================
% PLOT 1/V(omega,theta)
% ============================================================

Efun = matlabFunction( ...
    log(abs(1/Vtheta)), ...
    'Vars',[Om,theta]);

om_range = linspace(-4,4,200);
theta_range = linspace(-pi,pi,200);

[OM,TH] = meshgrid(om_range,theta_range);

E = Efun(OM,TH);

figure;

pcolor(OM,TH,E);
shading interp;

xlabel('\omega');
ylabel('\theta');

title('Energy plot: 1/V(\omega,\theta)');

colorbar;
% ============================================================
% FUNCTION 1:
% CANONICAL QUOTIENT BASIS
%
% Basis modulo c^2+s^2-1:
%
% {omega^i s^k}
% {omega^i c s^k}
%
% ============================================================

function Z = createQuotientBasis( ...
    n_omega,n_theta,omega,c,s)

    Z = [];

    for i = 0:2*n_omega

        % ----------------------------------------------------
        % Terms without c
        %
        % omega^i s^k
        % ----------------------------------------------------

        for k = 0:n_theta

            Z = [Z;
                 omega^i*s^k];

        end

        % ----------------------------------------------------
        % Terms containing one c
        %
        % omega^i c s^k
        %
        % c^2 is not needed because:
        %
        % c^2 = 1-s^2
        % ----------------------------------------------------

        for k = 0:n_theta-1

            Z = [Z;
                 omega^i*c*s^k];

        end

    end

end


% ============================================================
% FUNCTION 2:
% QUOTIENT-SPACE SOS GRAM BASIS
%
% This is the basis z used by sossosvar:
%
% s1v = z' Q z
%
% with z consisting of canonical quotient monomials.
% ============================================================

function Z = createQuotientSOSBasis( ...
    n_omega,n_theta,omega,c,s)

    Z = [];

    for i = 0:n_omega

        % omega^i s^k
        for k = 0:n_theta

            Z = [Z;
                 omega^i*s^k];

        end

        % omega^i c s^k
        for k = 0:n_theta-1

            Z = [Z;
                 omega^i*c*s^k];

        end

    end

end