function [m,m1] = match(u,v,l)
%% matching function & its partial derivative

% l is matching parameter (default = .407)
if (nargin < 3)  ||  isempty(l)
    l = optget('match','l',.407);   % matching parameter
end

m  = (u.*v)./(u.^l+v.^l).^(1/l);    % matching function
m1 = (m./u).^(1+l);                 % partial derivative

end