function [] = test_omp_thm2()
% test_OMP_Thm2 - Driver routine to test signal recovery algorithm.
%                 Specifically, test that OMP works using the value of N
%                 given in Theorem 2 of Tropp 2007      
%
% Syntax: 
%  [] = test_omp_thm2()
%
% Inputs:
%  None
%
% Outputs:
%  None
%
% Example:
%  test_omp_thm2()
%
% Dependencies:
%  plot_recovery
%  back_subs
% 
% TODO:
%
%
% Authors: JF,EY
% Revision history:
%  05 April 2014 - date written
%  08 April 2014 - I moved the subroutines to individual m-files.  
%                  test_recovery became plot_recovery.  The value of K
%                  can be bounded given the requirement delta>1/d (see
%                  "GC.pdf" page 7).
%  11 April 2014 - sparsity level m = ceil(0.95*d) is _very_ full; try
%                  m = ceil(0.05*d) for a sparse signal
                    


%% For repeatability, set PRNG (Mersenne twister) and seed (seed = 0)
rng('default');


%% General parameters
d = 50; % signal length
delta = 0.1; % 0 < delta < 0.36, 1-2*delta <= OMP recovery probability
K = 5; %For our particular choices, K<=6.7874 is good.  See GC.pdf


%% Generate reference signal and sparsify
m = ceil(0.05*d); % sparsity level (I think of m as amount filled)
% reference signal
s_full = 2*rand([d 1])-1; % uniform distribution on [-1,1]
num_remove_inds = d-m;
perm_inds = randperm(d);
remove_inds = perm_inds(1:num_remove_inds);
sparse_inds = perm_inds(num_remove_inds+1:d);
s = zeros([d 1]);
s(sparse_inds) = s_full(sparse_inds); % sparse reference signal


%% Measurment vectors
N = ceil(K*m*log(d/delta)); % N from Thm 2 of Tropp 2007

mu_Phi = zeros([N d]); % mean
mu_Sigma = 1/N*eye([d d]); % covariance
Phi = mvnrnd(mu_Phi,mu_Sigma); % measurement matrix, columns are vectors from
                               % multivariate N(0,1/N)
v = Phi*s; % data vector


%% Perform OMP

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
   [Q_Phi, R_Phi] = qrinsert(Q_Phi,R_Phi,t,Phi(:,lambda),'col');
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

% 7 - check to see if s_hat has
plot_recovery(s,s_hat)


end % test_omp_thm2

