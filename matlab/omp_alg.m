function [s_hat] = omp_alg(m,Phi,v)
% omp_alg - orthogonal matching pursuit algorithm
%
% Syntax: 
%  [s_hat] = tropp_fig1(m,Phi,v)
%
% Inputs:
%  m - sparsity level (integer)
%  Phi - data matrix (d by N matrix)
%  v - data vector (Phi*s, where s is original signal)
%
% Outputs:
%  s_hat - approximation to signal
%
% Dependencies:
%  back_subs
% 
% TODO:
%  Speed up the QR step
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written
                    

%% Perform OMP
[N,d] = size(Phi);

% 1 - initialization
resid = v; % residual vector
Lambda = []; % index set (column vector)
Phi_t = []; % matrix of atoms (will be filled by OMP)

Q_Phi = [];
R_Phi = [];

for t = 1:m
   
   % 2 - find index
   % if there are multiple indexes, max returns the first index
   [maxval,lambda] = max(abs(transpose(resid)*Phi));

   % 3 - augment index set and matrix of atoms
   Lambda = cat(1, Lambda, lambda); % add index column
   Phi_t = cat(2, Phi_t, Phi(:,lambda)); % add column on right

   % 4 - solve least squares problem
   % update QR decomp for new QR
   % TODO: qrinsert is implemented in pure matlab code
   %       it might be faster to just call QR (which is prolly 
   %       implemented in MKL)
   %[Q_Phi, R_Phi] = qrinsert(Q_Phi,R_Phi,t,Phi(:,lambda),'col');
   [Q_Phi,R_Phi] = qr(Phi_t);
   rhs = transpose(Q_Phi)*v;
   x = back_subs(R_Phi,rhs);
   
   %[q,r] = qr(Phi_t);
   %norm(Q_Phi-q,'fro')
   %norm(R_Phi-r,'fro')
   %norm(x-Phi_t\v)
  
   % 5 - update approximation and residual
   a = Phi_t*x;
   resid = v - a;
   
end

% 6 - construct approximate signal
s_hat = zeros([d 1]);
s_hat(Lambda) = x;

end % omp_alg
