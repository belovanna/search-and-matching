function [m,m1] = matchCD(u,v,l)
%% matching function in Cobb-Douglas form & its partial derivative

% l is matching parameter (default = .5)
if (nargin < 3)  ||  isempty(l)
    l = optget('match','l',.5);   % matching parameter
end

m  = (u.^l).*(v.^(1-l));    % matching function
m1 = l.*u.^(l-1).*v.^(1-l);  % partial derivative

end
