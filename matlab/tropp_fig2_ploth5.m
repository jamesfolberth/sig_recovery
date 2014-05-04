function [] = tropp_fig2_ploth5(h5file)
% tropp_fig2_ploth5 - routine to plot data (HDF5) for figure 1 on Tropp 2007
%
% Syntax: 
%  [] = tropp_fig2_ploth5()
%  [] = tropp_fig2_ploth5(h5file)
%
% Inputs:
%  h5file - path/to/fig2_data.h5
%
% Outputs:
%  None
%
% Examples:
%  fortran/tropp_fig2_comp  # run Fortran code
%  Using default save file
%  >> tropp_fig2_ploth5()
%
% Dependencies:
%  GNU Octave
%  fortran/*
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written
%  01 May 2014 - copied tropp_fig2_plot to tropp_fig2_ploth5
                    
% default savefile
if (nargin == 0)
   h5file = '../data/tropp_fig2_data_fort.h5';
end

% load data
load(h5file);

% plot it!
plot(m_vec,percent_recovered(:,:),'o-');
title(sprintf('Percentage of input signals recovered (d=%d) (Gaussian)',d));
xlabel('Sparsity level (m)');
ylabel('Percentage recovered');

legend(sprintf('N=%d',N_vec(1)),...
       sprintf('N=%d',N_vec(2)),...
       sprintf('N=%d',N_vec(3)),...
       sprintf('N=%d',N_vec(4)),...
       sprintf('N=%d',N_vec(5)),...
       'Location','NorthEast');

print('../figures/tropp_fig2.png','-dpng');

end % tropp_fig2_ploth5
