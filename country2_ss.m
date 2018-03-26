function [J,U,W,u,wp,t,p] = country2_ss(t,p,d,s,z)
%% Diamond-Mortensen-Pissarides (DMP) model for a two countries in a steady-state

% INPUT:
% t[2]    - steady-state response guessed values
% p[2]    - means of productivity
% sep[2]  - separation rates (default = [.0081,.0081])
% z[2]    - values of nonmarket activity (default = [.9,.8])
% d[2]    - discount of unemployment state abroad (default = [1.1,1.1])

% OUTPUT:
% J(p)     - values of filled job
% U(p)     - values of unemployment state
% W(p)     - values of employment state
% u(p)     - unemployment rates
% w(p)     - wage rates
% t(p)     - responses (market tightness)
% p        - states

%% 1. Set Parameters (Hagedorn and Manovskii, 2008, Table 2)

% globaly defined by optset
beta  = optget('country2','beta', 0.052);       % worker's bargaining power
delta = optget('country2','delta',.99^(1/12));  % discount rate
l     = optget('match','l',.407);               % matching parameter

% localy defined
if  ~exist('d', 'var') || isempty(d)           % discount of unemployment state abroad defaults
    d  = [1.1, 1.1];
end
if  ~exist('s', 'var') || isempty(s)           % separation rates defaults
    s  = [0.0081, 0.0081];
end
if  ~exist('z', 'var') || isempty(z)           % values of nonmarket activity defaults
    z  = [0.9, 0.8];
end

% Set the display format
%format long   % long fixed decimal
format shortG % short fixed decimal or scientific, with 5 digits

%% 2. compute certainty-equivalent steady-state
% [theta p] guessed values.
% NOTE : broyden computes Jacobian around these point.
% ---- : values should be close enough to solution, otherwise broyden will
% go out of positive values for t1,t2 during the search.
x0 = [t(1), p(1), t(2), p(2)]';
% the means of productivities
p0 = p';

f = @(x) DMP2_ss(x,p0,s,z,d);
optset('broyden','showiters','true');
%optset('broyden','tol',1e-16);
%optset('broyden','maxit',50);
[xstar,fval,flag] = broyden(f,x0); % Broyden's Inverse Method

if flag ~= 0
    warning('Broyden:FailToConverge. fval = "[%f %f %f %f]" ', fval);
end

%% 3. steady-state values
% responses
t = [xstar(1), xstar(3)];
% productivities (states)
p = [xstar(2), xstar(4)];

% cost of vacancy posting
c = cost(p);
% probability vacancy being filled next period
q = probvacancy(t,l);

% unemployment states cases ------------------------------
U0 = z+beta/(1-beta)*c.*t;
D0 = [U0(2)/d(1)-U0(1), U0(1)/d(2)-U0(2)];

D = D0/(1-delta);
for i=0:100
    D_tmp = D;
    
    h = heaviside(D); % step function
    
    D(1) = D0(1)/(1-delta*(1-h(1))) + delta/(1-delta)*D0(2)/d(1)*h(2);
    D(2) = D0(2)/(1-delta*(1-h(2))) + delta/(1-delta)*D0(1)/d(2)*h(1);
    
    dif = sum( (D-D_tmp).^2 );
    %disp([i,h1,h2,D1,D2,dif]);
    
    if dif < eps
        break
    end
end

% unemployment states
U = U0/(1-delta) + delta/(1-delta)*D.*h;

% check case 1: U1>U2/d1 (D1<0), U2>U1/d2 (D2<0) - no migration
mcase = 1;
if ( D(1)>0 || D(2)>0 )
    % check case 2 : U1>U2/d1 (D1<0), U2<U1/d2 (D2>0) - 2=>1
    mcase = 2;
    if ( D(2)<0 )
        % check case 3 : U1<U2/d1 (D1>0), U2>U1/d2 (D2<0) - 1=>2
        mcase = 3;
    end
    if ( D(1)>0 && D(2)>0 )
        error('Unemployment error! No mcase found!');
    end
end

fprintf('\nmcase = %i: D(1),D(2) = %f, %f.\n',mcase,D);

% wages
wp = p + (delta*(1-s)-1).*c./(delta*q);

% unemployment rates
u = s./(s+t.*q);

% filled job values
J = c./(delta*q);

% employment state values
W = beta/(1-beta)*J + U + D.*h;

end