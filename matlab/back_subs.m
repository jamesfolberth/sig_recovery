function [x] = back_subs(R,b)
% back_subs - perform backward substitution on an upper trapezoidal matrix
%             this algorithm assumes that R is ''tall and skinny'', meaning
%             that num_rows >= num_cols.
%
% Syntax: 
%  [x] = back_subs(R,b)
%
% Inputs:
%  R - upper trapezoidal (triangular) matrix (e.g. from QR decomp)
%  b - RHS of linear system R*x=b
%
% Outputs:
%  x - solution of system
%
% Example:
%  A = rand([10 5]);
%  rhs = rand([10 1]);
%  [Q,R] = qr(A);
%  b = transpose(Q)*rhs;
%  x = back_subs(R,b);
%
% Dependencies:
%  None

ZERO_TOL = 10^(-14);

[m,n] = size(R);
minmn = min([m,n]);

x = zeros([n 1]);

if (abs(R(minmn,minmn)) < ZERO_TOL)
   error('back_subs: input matrix R is singular to precision ZERO_TOL');
end

x(minmn) = b(minmn)/R(minmn,minmn);

for i=minmn-1:-1:1
   if (abs(R(i,i)) < ZERO_TOL)
      error('back_subs: input matrix R is singular to precision ZERO_TOL');
   end

   x(i) = (b(i) - R(i,i+1:minmn)*x(i+1:minmn))/R(i,i);
end

end % back_subs
