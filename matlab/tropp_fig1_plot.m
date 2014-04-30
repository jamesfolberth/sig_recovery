function [] = tropp_fig1_plot(matfile)
% tropp_fig1_plot - routine to plot data for figure 1 on Tropp 2007
%
% Syntax: 
%  [] = tropp_fig1_plot()
%  [] = tropp_fig1_plot(matfile)
%
% Inputs:
%  matfile - path/to/fig1_data.mat
%
% Outputs:
%  None
%
% Examples:
%  Using default save file
%  >> tropp_fig1_comp()
%  >> tropp_fig1_plot()
%
%  Using given save file
%  >> tropp_fig1_comp('savefile.mat')
%  >> tropp_fig1_plot('savefile.mat')
%
% Dependencies:
%  omp_alg
%  check_recovery
%  tropp_fig1_comp
%
% TODO:
%  Test if default save file works for mac
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written
                    
% default savefile
if (nargin == 0)
   if strcmp(computer(),'GLNXA64')
      matfile = '../data/tropp_fig1_data_glnx64.mat';
   elseif strcmp(computer(),'MACI64')
      matfile = '../data/tropp_fig1_data_mac64.mat';
   else
      error('tropp_fig1_plot: default file not found!');
   end
end

% load data
load(matfile, '-mat');

% plot it!
plot(N_vec,percent_recovered(:,:),'o-');
title(sprintf('Percentage of input signals recovered (d=%d) (Gaussian)',d));
xlabel('Number of measurments (N)');
ylabel('Percentage recovered');

legend(sprintf('m=%d',m_vec(1)),...
       sprintf('m=%d',m_vec(2)),...
       sprintf('m=%d',m_vec(3)),...
       sprintf('m=%d',m_vec(4)),...
       sprintf('m=%d',m_vec(5)),...
       'Location','SouthEast');

% TODO call matlab2tikz or whatever

end % tropp_fig1_plot
