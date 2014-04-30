function [] = plot_recovery(x,y)
% plot_recovery - plot the reference signal, the approximation,
%                 and their absolute difference
%
% Syntax: 
%  [] = plot_recovery(s,v)
%
% Inputs:
%  x - reference signal; d vector
%  y - sparsified signal; d vector
%
% Outputs:
%  None
%
% Example:
%  plot_recovery([1;2;3;4],[1.01;1.97;3.02;4.01]);
%
% Dependencies:
%  None

subplot(3,1,1)
plot(x)
%ylabel('ref sig')

subplot(3,1,2);
plot(y)
%ylabel('approx sig')

subplot(3,1,3);
plot(abs(x-y))
%ylabel('abs diff')
%xlabel('index')

end % plot_recovery
