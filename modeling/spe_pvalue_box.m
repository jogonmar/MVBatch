function p_value = spe_pvalue_box(res,lim)

% Inverse of the control limit for SPE or Q statistic based on P Nomikos,
% JF MacGregor. Multivariate SPC charts for monitoring batch processes. Technometrics 37 (1), 41-59.
%
% p_value = spe_pvalue(res,lim,pcs_left)       % complete call
%
% INPUTS:
%
% res: (NxM) Two-way residuals data matrix, N(observations) x M(variables)
%
% lim: (1x1) control limit at a 1-p_value confidence level.
%
%
% OUTPUTS:
%
% p_value: (1x1) p-value of the test.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: 10/Apr/16
%
% Copyright (C) 2016  José M. González Martínez
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
if ndims(res)~=2, error('Incorrect number of dimensions of res.'); end;
s = size(res);
if find(s<1), error('Incorrect content of res.'); end;
if (lim<0), error('Incorrect value of lim.'); end;

% Main code
m = mean(sum(res.^2,2));
v = var(sum(res.^2,2));
p_value = 1-chi2cdf(lim/(v/(2*m)),(2*m^2)/v);
 


