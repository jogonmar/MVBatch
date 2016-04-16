function n_pa=number_params_cov(phases,n_var)

% Number of parameters of the covariance-matrix.
%
% n_pa=number_params_cov(phases,n_var) % complete call
%
% INPUTS:
%
% phases: (n_phasesx5) original set of disjoints phases modelling an interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
%
% n_var: (1x1) number of process variables.
%
%
% OUTPUTS:
%
% n_pa: (1x1) number of parameters. 
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 19/Abr/07.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;

% Main code

sp=size(phases);

n_pa=0;
for i=1:sp(1),   
  var=n_var*(phases(i,3)+1); 
  mod = ((var+1)*var)/2;  
  n_pa = n_pa + mod*phases(i,2);
end
    