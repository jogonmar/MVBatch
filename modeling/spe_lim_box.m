function lim = spe_lim_box(res,p_value)

% Control limit for SPE or Q statistic based on the Box approximation published in P Nomikos,
% JF MacGregor. Multivariate SPC charts for monitoring batch processes. Technometrics 37 (1), 41-59.
%
% lim = spe_lim_box(res,p_value)        % complete call
%
% INPUTS:
%
% res: [NxM] Two-way residuals data matrix
%
% p_value: [1x1] p-value of the test, in (0,1]
%
%
% OUTPUTS:
%
% lim: [1x1] control limit at a 1-p_value confidence level.
%
%
% coded by: José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2017  José M. González Martínez
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

% Validate dimensions of input data
assert (isequal(size(p_value), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (p_value>=0 && p_value<1, 'Value Error: 2nd argument must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);


%% Main code
m = mean(sum(res.^2,2));
v = var(sum(res.^2,2));
lim = (v/(2*m))*chi2inv(1-p_value,(2*m^2)/v);

