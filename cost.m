%function [c,c1] = cost(cK,cW,xi)
function [c,c1] = cost(xi)
%% cost of vacancy posting & its first derivative
xi=[0.2,0.2];
% default arguments
if ~exist('cK','var') || isempty(cK)
    cK = optget('cost','cK',.474);
end
if ~exist('cW','var') || isempty(cW)
    cW = optget('cost','cW',.110);
end
if ~exist('xi','var') || isempty(xi)
    xi = optget('cost','xi',.449);
end
%c  = cK.*p + cW.*p.^xi;
%c1 = cK + cW.*xi.*p.^(xi-1);
c=xi;
c1=[0,0];
end