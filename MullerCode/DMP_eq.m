function [out1,out2,out3] = DMP_eq(flag,s,x,Eh,e,beta,delta,rho,sep,z,c,c1,q,q1)

[ns dx] = size(x);          % number of nodes, dimension of response vector
ds      = size(s,2);        % dimension of state vector
dh      = size(Eh,2);       % dimension of expectational vector

% responses
theta = x;                  % theta denotes log(theta)  

% states
p = s(:,1);                 % denotes log(p)


switch flag
    case 'b';
        % lower bound
        out1      = -inf*ones(ns,dx);
        
        % upper bound
        out2      = inf*ones(ns,dx);
        
    case 'f';
        % equilibrium function
        out1      = zeros(ns,dx);
        out1(:,1) = Eh-c(exp(p))./(delta*q(exp(theta)));
        
        % derivative with respect to x
        out2        = zeros(ns,dx,dx);
        out2(:,1,1) = c(exp(p)).*q1(exp(theta)).*(exp(theta))./(delta*q(exp(theta)).^2);    %theta
        
        %derivative with respect to Eh
        out3        = zeros(ns,dx,dh);
        out3(:,1,1) = ones(ns,1);               %E[.]
        
    case 'g';
        % transition equation
        out1      = zeros(ns,ds);
        out1(:,1) = rho*p+e;
        
        % derivative with respect to x
        out2      = zeros(ns,ds,dx);
        
    case 'h';
        % expectation function
        out1      = zeros(ns,dh);
        out1(:,1) = (1-beta)*(exp(p)-z)-c(exp(p)).*(beta*exp(theta))+(1-sep)*c(exp(p))./q(exp(theta));
        
        % derivative with respect to x
        out2        = zeros(ns,dh,dx);
        out2(:,1,1) = -c(exp(p)).*(beta*exp(theta))-(1-sep)*c(exp(p)).*q1(exp(theta)).*(exp(theta))./q(exp(theta)).^2;        %theta
        
        % derivative with respect to s
        out3        = zeros(ns,dh,ds);
        out3(:,1,1) = (1-beta)*exp(p)-c1(exp(p)).*(beta*exp(theta)).*exp(p)+(1-sep)*c1(exp(p)).*(exp(p))./q(exp(theta));    %p
end


end

