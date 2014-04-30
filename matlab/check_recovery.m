function [recovered,nrm] = check_recovery(s,s_hat)
% check_recovery - determine if s_hat (the approximation) matches
%                  the original signal s
%
% Syntax: 
%  [recovered] = check_recovery(s,s_hat)
%  [recovered, nrm] = check_recovery(s,s_hat)
%
% Inputs:
%  s - original signal (d vector)
%  s_hat - approximate signal (d vector)
%
% Outputs:
%  recovered - 1 if recovered, 0 if not
%  nrm - norm of difference (optional)
%
% Dependencies:
%  None
% 
% TODO:
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written

d = numel(s);
ZERO_TOL = 10^(-14);
RECOVERY_TOL = ZERO_TOL;

recovered = 0;
nrm_ = norm(s-s_hat,inf);
if (nrm_ < RECOVERY_TOL)
   recovered = 1;
end

if (nargout > 1)
   nrm = nrm_;
end

end % check_recovery
