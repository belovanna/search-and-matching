  function y = DMP2_ss(x,p0,s,z,d)
%% Diamond-Mortensen-Pissarides (DMP) model at certainty-equivalent steady-state

%disp([x(1),x(3)]); % iteration's input. NOTE: they have to be positive for convergence

% Default Parameters
beta  = optget('country','beta', 0.052);       % worker's bargaining power
delta = optget('country','delta',.99^(1/12));  % discount rate
l     = optget('match','l',.407);              % matching parameter

% responses
t = [x(1), x(3)];      % theta denotes theta(p)!
% states
p = [x(2), x(4)];      % p denotes p!

% cost of vacancy posting
c = cost(p);
% probability vacancy being filled next period. NOTE: t should be positive!
q = probvacancy(t,l);

%% Unemployment States
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

%D = 0; % fully independent solution

%% Steady State Conditions

% rational expectations equation
y(1,1) = (1-beta)*(p(1)-z(1))...
    - beta*c(1)*t(1) + (1-s(1))*c(1)/q(1) - c(1)/(delta*q(1))...
    - (1-beta)*D(1)*h(1);
y(3,1) = (1-beta)*(p(2)-z(2))...
    - beta*c(2)*t(2) + (1-s(2))*c(2)/q(2) - c(2)/(delta*q(2))...
    - (1-beta)*D(2)*h(2);
% productivity
y(2,1) = p(1)-p0(1);
y(4,1) = p(2)-p0(2);
%{
% unemployment rate
y(3,1) = (1-delta)*U(1)...
    - z(1) - beta/(1-beta)*c(1)*t(1)...
    - delta*D(1)*h(1);
y(6,1) = (1-delta)*U(2)...
    - z(2) - beta/(1-beta)*c(2)*t(2)...
    - delta*D(2)*h(2);
%}

end