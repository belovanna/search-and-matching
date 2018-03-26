function y = DMP_ss(x,p0,s,z)
%% Diamond-Mortensen-Pissarides (DMP) model at certainty-equivalent steady-state

%% 1. set parameters
% globally defined parameters
beta  = optget('country','beta', 0.052);        % worker's bargaining power
delta = optget('country','delta', .99^(1/12));  % discount rate
l     = optget('match','l', .407);              % matching parameter
% localy defined parameters
if  ~exist('s', 'var') || isempty(s)
    s = optget('country','s', .0081);           % separation rate
end
if  ~exist('z', 'var') || isempty(z)
    z = optget('country','z', .955);            % value of nonmarket activity
end

%% 2. responses & states
% responses
theta = x(1,1);        % theta denotes theta(p)!
% states
p     = x(2,1);        % p denotes p!

c = cost(p);                % cost of vacancy posting
q = probvacancy(theta,l);   % probability vacancy being filled next period

%% 3. steady state conditions
% rational expectations equation
y(1,1) = (1-beta)*(p-z) - c*beta*theta + (1-s)*c/q - c/(delta*q);
% productivity
y(2,1) = p-p0;

end