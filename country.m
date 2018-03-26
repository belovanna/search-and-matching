function [J,U,W,u,wp,theta,p,model] = country(nshocks,rho,sig);
%% Diamond-Mortensen-Pissarides (DMP) model for a single country
% INPUT:
% nshocks - number of discretized nodes (default = 5)
% rho     - persistence of productivity process (default = .9895)
% sig     - variance of innovations in productivity process (default = .0034)
% OUTPUT:
% J(p)     - value of filled job
% U(p)     - value of unemployment state
% W(p)     - value of employment state
% u(p)     - unemployment rate
% w(p)     - wage rate
% theta(p) - responce (market tightness)
% p        - state
% model    - model structure to be used in DMP_simulation

%% 1. Set country parameters (Hagedorn and Manovskii, 2008, Table 2)

if ~exist('nshocks','var') || isempty(nshocks)
    nshocks = 5;           % number of nodes
end
if ~exist('rho','var') || isempty(rho)
    rho = .9895;           % persistence of productivity process
    optset('country','rho',rho);
end
if ~exist('sig','var') || isempty(sig)
    sig = .0034;           % variance of innovations in productivity process
end

% function parameters (globaly defined by optset)
beta  = optget('country','beta', 0.052);       % worker's bargaining power
delta = optget('country','delta',.99^(1/12));  % discount rate
sep   = optget('country','sep',  .0081);       % separation rate
z     = optget('country','z',    .955);        % value of nonmarket activity

%% 2. Compute certainty-equivalent steady-state
% guessed values
theta0 = 0.6;  
p0     = 1.0;
% [job_value,unemployment_value,employment_value,unemployment_rate,wage,responce,state]
[J_ss,U_ss,W_ss,u_ss,wp_ss,theta_ss,p_ss] = country_ss(theta0,p0);

%% 3. Discretization of the innovation in the productivity process
[e,w] = qnwnorm(nshocks,0,sig.^2);      % normal innovations

%% 4. Set up model structure
clear model;
%model.func = 'DMP_mfunc';  % model function file
model.func = 'DMP_mfuncU'; % model function file
model.e = e;               % shocks
model.w = w;               % corresponded probabilities
model.params = {U_ss};         % we passed parameters using global optset

%% 5. Discretization of the state space, log(p)
n    = 4; % number of nodes
smin = log(p_ss)+min(e);
smax = log(p_ss)+max(e);

% Defines a function family structure used to define state space nodes (1-Dimensional).
% ---- Chebychev collocation.
% 0 - Gaussian nodes (do not include endpoints) - default
% 1 - Gaussian nodes extended to endpoints
% 2 - Lobatto nodes  (includes endpoints)
optset('chebnode','nodetype',0);
fspace = fundefn('cheb',n,smin,smax);
% ---- Spline collocation (includes endpoints).
%fspace = fundefn('spli',n,smin,smax,3);              % 3 - cubic splines (default)

% Computes default nodes for a family of functions (n-vector)
scoord = funnode(fspace);
% Expands matrices into the associated grid points (cell array)
snodes = gridmake(scoord); % NOTE: due to 1D-state problem, snodes=scoord.

%% 6. Initial guess (CE steady-state)
c = cost(p_ss);            % cost of vacancy posting
q = probvacancy(theta_ss); % probability vacancy being filled
% initial guess for equilibrium responses at nodes: log(theta),U.
xinit = repmat([log(theta_ss),1],size(snodes,1),1);
% initial guess for expectation variable at nodes: E[.].
hinit = repmat([(1-beta)*(p_ss-z)-c*beta*theta_ss+(1-sep)*c/q,delta],size(snodes,1),1);

%% 7. Solve the model

% --- OPTIONS FOR REMSOLVE (rational expectations model solver)
optset('remsolve','maxit',10000);      % maximum number of iterations (500 is default)
%optset('remsolve','tol',10^-12);      % convergence tolerance (sqrt(eps) is default)
%optset('remsolve','nres',10);         % nres*fspace.n uniform nodes to evaluate residual (default is 10)
optset('remsolve','showiters',1);     % 0/1, 1 to display iteration results (defaults)
% --- OPTIONS FOR ARBIT (arbitrage equation solver at specified nodes)
%optset('arbit','maxit',1000);         % maximum number of iterations for CP (300 is default)
%optset('arbit','tol',10^-12);         % convergence tolerance for CP (sqrt(eps) is default)
%optset('arbit','lcpmethod','smooth'); % 'minmax' (default) or 'smooth'

% OUTPUT:
% coeff - expectation function approximation basis coefficients
% s     - residual evaluation coordinates
% x     - equilibrium responses at evaluation points
% h     - expectation variables at evaluation points
% f     - marginal arbitrage benefits at evaluation points
% resid - rational expectation residuals at evaluation points
%[coeff,s,x,h,f,resid] = remsolve(model,fspace,scoord,xinit,hinit);
[coeff,s,x,h,f,resid] = remsolve(model,fspace,scoord,xinit,hinit);

% states            
p = exp(s);
% responses
theta = exp(x(:,1));
U     = x(:,2)*U_ss;

c     = cost(p);                            % cost of vacancy posting
q     = probvacancy(theta);                 % probability vacancy being filled next period

u     = sep./(sep+match(1,theta));          % unemployment rate
wp    = beta*p+(1-beta)*z+c.*(beta*theta);  % wage rate

J     = p-wp+delta*(1-sep)*c./(delta*q);
WU    = beta/(1-beta) * J;                  % W-U
W = U + WU;

%% 8. plots

f=figure(1);
hold on;

% Plot Chebychev residual (surface)
subplot(4,2,1);
h_resid=plot(p,resid);
%title('Chebychev Residual');
xlabel('p, productivity');
ylabel('Chebychev Residual');
xlim([min(p) max(p)]);

% Plot responses
subplot(4,2,2);
h_theta=plot(p,theta);
xlabel('p, productivity');
ylabel('theta, market tightness');
xlim([min(p) max(p)]);

% Plot wage rates
subplot(4,2,3);
h_wp=plot(p,wp);
xlabel('p, productivity');
ylabel('w(p), wage rate');
xlim([min(p) max(p)]);

% Plot unemployment rate
subplot(4,2,4);
h_u=plot(p,u);
xlabel('p, productivity');
ylabel('u, unemployment rate');
xlim([min(p) max(p)]);

subplot(4,2,5);
h_U=plot(p,U);
xlabel('p, productivity');
ylabel('U, unemployment flow');
xlim([min(p) max(p)]);

subplot(4,2,6);
h_W=plot(p,W);
xlabel('p, productivity');
ylabel('W, employment flow');
xlim([min(p) max(p)]);

subplot(4,2,7);
h_WU=plot(p,WU);
xlabel('p, productivity');
ylabel('W-U flow');
xlim([min(p) max(p)]);

subplot(4,2,8);
h_J=plot(p,J);
xlabel('p, productivity');
ylabel('J, filled jobs flow');
xlim([min(p) max(p)]);

hold off;

end
