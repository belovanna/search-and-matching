function [out1,out2,out3] = DMP_mfunc(flag,s,x,Eh,e)
%% Diamond-Mortensen-Pissarides (DMP) model function
% Non-Linear Rational Expectations function, ready for remsolve.

% INPUT:
%
% flag='b' - returns bound function
%   xl:n.dx, xu:n.dx
% flag='f' - returns reward function and derivatives
%   f:n.dx, fx:n.dx.dx, feh:n.dx.dh
% flag='g' - returns transition function and derivatives
%   g:n.ds, gx:n.ds.dx
% flag='h' - returns expectation function and derivatives
%   h:n.dh, hx:n.dh.dx, hs:n.dh.ds
%
% s - states, log(p)
% x - initial guess for equilibrium responses at nodes, log(theta)
% Eh - expectation vector
% e - shocks

% function parameters (globaly defined by optset)
beta  = optget('country','beta', 0.052);       % worker's bargaining power
delta = optget('country','delta',.99^(1/12));  % discount rate
rho   = optget('country','rho',.9895);         % persistence of productivity process
sep   = optget('country','sep',  .0081);       % separation rate
z     = optget('country','z',    .955);        % value of nonmarket activity
l     = optget('match','l',.407);              % matching parameter

[n,dx] = size(x);          % number of collocation nodes, dimension of response vector
ds     = size(s,2);        % dimension of state vector
dh     = size(Eh,2);       % dimension of expectational vector

% responses
theta = exp(x); 
% states
p = exp(s(:,1));

[c,c1] = cost(p);              % cost of vacancy posting
[q,q1] = probvacancy(theta,l); % probability vacancy being filled next period

switch flag
    case 'b';
        % lower bound
        out1      = -inf*ones(n,dx);
        
        % upper bound
        out2      = inf*ones(n,dx);
        
    case 'f';
        % equilibrium function
        out1      = zeros(n,dx);
        out1(:,1) = Eh-c./(delta*q);
        
        % derivative with respect to x
        out2        = zeros(n,dx,dx);
        out2(:,1,1) = theta.*c.*q1./(delta*q.^2);    % ln(theta)
        
        %derivative with respect to Eh
        out3        = zeros(n,dx,dh);
        out3(:,1,1) = ones(n,1);                     % E[.]
        
    case 'g';
        % transition equation
        out1      = zeros(n,ds);
        out1(:,1) = rho*log(p)+e;
        
        % derivative with respect to x
        out2      = zeros(n,ds,dx);                  % ln(theta)
        
    case 'h';
        % expectation function
        out1      = zeros(n,dh);
        out1(:,1) = (1-beta)*(p-z)-c.*(beta*theta)+(1-sep)*c./q;
        
        % derivative with respect to x
        out2        = zeros(n,dh,dx);
        out2(:,1,1) = -theta.*(c*beta+(1-sep)*c.*q1./q.^2);          % ln(theta)
        
        % derivative with respect to s
        out3        = zeros(n,dh,ds);
        out3(:,1,1) = p.*((1-beta)-c1.*(beta*theta)+(1-sep)*c1./q);  % ln(p)
end

end

