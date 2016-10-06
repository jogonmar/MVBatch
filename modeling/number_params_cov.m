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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Abr/07
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    