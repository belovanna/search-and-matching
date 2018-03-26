% Replicates Hagedorn and Manovskii (AER, 2008).

%% 0. housekeeping
clear all;
clc;

%% 1. set parameters (Hagedorn and Manovskii, 2008, Table 2)

% parameters
beta  = 0.052;           % worker's bargaining power
delta = .99^(1/12);      % discount rate
l     = .407;            % matching parameter
rho   = .9895;           % persistence of productivity process
sep   = .0081;           % separation rate
sig   = .0034;           % variance of innovations in productivity process
z     = .955;            % value of nonmarket activity

% functional forms
m  = @(u,v) (u.*v)./(u.^l+v.^l).^(1/l);         % matching function
m1 = @(u,v) (m(u,v)./u).^(1+l);                 % partial derivative

c  = @(p) .474*p+.110*p.^(.449);                % cost of vacancy posting
c1 = @(p) .474+.110*.449*p.^(.449-1);

q  = @(theta) m(1./theta,1);                    % probability vacancy being filled next period
q1 = @(theta) m1(1./theta,1).*(-1./theta.^2);

%% 2. compute certainty-equivalent steady-state
% [log(theta) log(p)]
x0 = [log(.64) log(1)]';
%optset('broyden','showiters','true');
f = @(x) DMP_ss(x,beta,delta,sep,z,c,q);
[xstar,fval] = broyden(f,x0);

% steady-state values
thetastar = xstar(1,1);                     % theta* denotes log(theta*)!
pstar     = xstar(2,1);                     % p* denotes log(p*)!
ustar     = sep/(sep+m(1,exp(thetastar)));

%% 3. discretization of the innovation in the productivity process
nshocks = 5;                            % number of nodes
[e,w] = qnwnorm(nshocks,0,sig.^2);      % normal innovations

%% 4. set up model structure
clear model
model.func = 'DMP_eq';
model.e = e;
model.w = w;
model.params = {beta,delta,rho,sep,z,c,c1,q,q1};

%% 5. discretization of the state space, log p
n     = 4;
smin = pstar+min(e);
smax = pstar+max(e);

fspace = fundefn('cheb',n,smin,smax);               % Chebychev collocation, could also use splines, 'spli'.
scoord = funnode(fspace);
snodes = gridmake(scoord);

%% 6. initial guess (CE steady-state)
% response, theta.
xinit = repmat(thetastar,size(snodes,1),1);
% expectations, E[.].
hinit = repmat((1-beta)*(exp(pstar)-z)-c(exp(pstar))*beta*exp(thetastar)+(1-sep)*c(exp(pstar))/q(exp(thetastar)),size(snodes,1),1);

%% 7. solve the model
optset('remsolve','maxit',1000);
%optset('remsolve','tol',10^-12);
[coeff,s,x,h,f,resid] = remsolve(model,fspace,scoord,xinit,hinit);

% states            
p = s;                   % p denotes log(p)

% response
theta = x;               % theta denotes log(theta)

%% 8. simulation 

nweeks = 3600;           
npath  = 200;
sinit  = pstar*ones(npath,1);

[spath,xpath] = remsimul(model,sinit,nweeks,s,x);

% log(theta)
thetapath = xpath;
% log(p)
ppath     = spath;
% log(u)
upath = zeros(npath,nweeks+1);
upath(:,1) = ustar;
for t=2:nweeks+1
    upath(:,t)=1-((1-sep)*(1-upath(:,t-1))+m(1,exp(thetapath(:,t-1))).*upath(:,t-1));
end
upath = log(upath);
% log(v)
vpath = upath+thetapath;

%% 9. select and HP-filter time paths

upath_q = upath(:,1:12:nweeks+1);                   % pick quarterly data
upath_f = hpfilter(upath_q,1600)';                  % HP-filter quarterly data
upath_d = upath_q-upath_f;                          % log-deviation from HP-trend

vpath_q = vpath(:,1:12:nweeks+1);
vpath_f = hpfilter(vpath_q,1600)';
vpath_d = vpath_q-vpath_f;  

thetapath_q = thetapath(:,1:12:nweeks+1);
thetapath_f = hpfilter(thetapath_q,1600)';
thetapath_d = thetapath_q-thetapath_f;  

ppath_q = ppath(:,1:12:nweeks+1);
ppath_f = hpfilter(ppath_q,1600)';
ppath_d = ppath_q-ppath_f;  

%% 10. replicate Hagedorn and Manovskii (2008, Table 4)

% average standard deviation
ustd     = mean(std(upath_d,0,2));
vstd     = mean(std(vpath_d,0,2));
thetastd = mean(std(thetapath_d,0,2));
pstd     = mean(std(ppath_d,0,2));

fprintf('\n\t u \t \t \t v \t \t v/u \t \t p\n');
disp('Standard deviation:')
disp([ustd vstd thetastd pstd]);

% average autocorrelation
auto = zeros(4,size(upath_q,1));
for i=1:size(auto,2)-1
    corrmat   = corrcoef([upath_d(i,1:end-1)',upath_d(i,2:end)']);
    auto(1,i) = corrmat(1,2);
    corrmat   = corrcoef([vpath_d(i,1:end-1)',vpath_d(i,2:end)']);
    auto(2,i) = corrmat(1,2);
    corrmat   = corrcoef([thetapath_d(i,1:end-1)',thetapath_d(i,2:end)']);
    auto(3,i) = corrmat(1,2);
    corrmat   = corrcoef([ppath_d(i,1:end-1)',ppath_d(i,2:end)']);
    auto(4,i) = corrmat(1,2);
end

auto_mat = mean(auto,2)';
disp('Quarterly autocorrelation:')
disp(auto_mat);

% average correlation
corr = zeros(4,4,size(upath_q,1));
for i=1:size(corr,3)
    corr(:,:,i)= corrcoef([upath_d(i,:)',vpath_d(i,:)',thetapath_d(i,:)',ppath_d(i,:)']);
end

cor_mat = triu(mean(corr,3));
disp('Correlation matrix:')
disp(cor_mat);

%% 11. plot residuals

% plot Chebychev residual (surface)
 figure(1)
 hh=plot(exp(s),resid);
 %title('Chebychev Residual')
 xlabel('p, productivity'); 
 ylabel('Chebychev Residual');