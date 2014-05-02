function [] = tropp_fig1_ploth5(h5file)
% tropp_fig1_ploth5 - routine to plot data (HDF5) for figure 1 on Tropp 2007
%
% Syntax: 
%  [] = tropp_fig1_ploth5()
%  [] = tropp_fig1_ploth5(h5file)
%
% Inputs:
%  h5file - path/to/fig1_data.h5
%
% Outputs:
%  None
%
% Examples:
%  fortran/tropp_fig1_comp  # run Fortran code
%  Using default save file
%  >> tropp_fig1_ploth5()
%
% Dependencies:
%  GNU Octave
%  fortran/*
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written
%  01 May 2014 - copied tropp_fig1_plot to tropp_fig1_ploth5
                    
% default savefile
if (nargin == 0)
   h5file = '../data/tropp_fig1_data_fort.h5';
end

% load data
load(h5file);

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

end % tropp_fig1_ploth5
