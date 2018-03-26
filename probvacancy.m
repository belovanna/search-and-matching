function [q,q1] = probvacancy(theta,l)
%% probability vacancy being filled next period

if theta <= 0
    error('Probvacancy error! Theta = %d is negative (or zero)!', theta);
end

% l is matching parameter (default = .407)
if (nargin < 2)  ||  isempty(l)
    l = optget('matchfunc','l',.407);        % matching parameter
end

[q,q1] = match(1./theta,1,l);
q1 = q1.*(-1./theta.^2);

end

