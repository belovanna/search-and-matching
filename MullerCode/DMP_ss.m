function y = DMP_ss(x,beta,delta,sep,z,c,q)
% responses
theta = x(1,1);                 % theta denotes theta(p)!

% states
p = x(2,1);                     % p denotes log(p)!

% steady state conditions

% rational expectations equation
y(1,1) = (1-beta)*(exp(p)-z)-c(exp(p))*beta*exp(theta)+(1-sep)*c(exp(p))/q(exp(theta))-c(exp(p))/(delta*q(exp(theta)));
% productivity
y(2,1) = exp(p)-1;

end