% Function to compute Eccentric anomaly from Mean anomaly using Kepler
% equation. Iteration will stop when difference between two consecutive
% iteration are smaller than tol value.
%
% Input:  e - numerical eccentricity of ellipse (---),
%         Mk - mean anomaly in (rad),
%         tol - tolerance to stop iterative process in arcseconds (").
%
% Output: Ek - eccentric anomaly in (rad).
%
%% ------------------------------------------------------------------------

function Ek = kepler(e, Mk, tol)
E = [Mk, Mk + e*sin(Mk)];

i = 2;
while abs(E(i) - E(i-1)) > tol*(pi/648000) 
    i = i + 1;
    E(i) = Mk + e*sin(E(i-1));
end

Ek = E(end);
