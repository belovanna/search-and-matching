function [J,U,W,u,wp,theta,p] = country_ss(theta,p,s,z)
%% Diamond-Mortensen-Pissarides (DMP) model for a single country in a steady-state
% INPUT:
% theta - response guessed value
% p     - state guessed value
% s     - separation rate
% z     - value of nonmarket activity
% OUTPUT:
% J(p)     - value of filled job
% U(p)     - value of unemployment state
% W(p)     - value of employment state
% u(p)     - unemployment rate
% w(p)     - wage rate
% theta(p) - response (market tightness)
% p        - state

%% 1. set parameters (Hagedorn and Manovskii, 2008, Table 2)
% globaly defined by optset
beta  = optget('country','beta', 0.052);        % worker's bargaining power
delta = optget('country','delta', .99^(1/12));  % discount rate
l     = optget('match','l', .407);              % matching parameter
% localy defined
if  ~exist('s', 'var') || isempty(s)
    s = optget('country','s', .0081);           % separation rate default
end
if  ~exist('z', 'var') || isempty(z)
    z = optget('country','z', .955);            % value of nonmarket activity default
end

%% 2. compute certainty-equivalent steady-state
% [theta p] guessed values
x0 = [theta p]';
% the mean of productivity (for single country usually normalized to one)
if  ~exist('p', 'var') || isempty(p)
    p0 = 1;
else
    p0 = p;
end

f = @(x) DMP_ss(x,p0,s,z);
%optset('broyden','showiters','true');
[xstar,fval,flag] = broyden(f,x0); % Broyden's Inverse Method

if flag ~= 0
    warning('Broyden:FailToConverge. fval = "[%f %f]" ', fval);
end

%% 3. steady-state values
theta = xstar(1,1);                     % response
p     = xstar(2,1);                     % state

c     = cost(p);                        % cost of vacancy posting
q     = probvacancy(theta,l);           % probability vacancy being filled next period

u     = s/(s+match(1,theta));           % unemployment rate
wp    = beta*p+(1-beta)*z+c*beta*theta; % wage rate

J     = c/(delta*q);                            % filled job value
W     = (wp-delta*beta/(1-beta)*s*J)/(1-delta); % employment state value
U     = (z+c*beta/(1-beta)*theta)/(1-delta);    % unemployment state value

end
