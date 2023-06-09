function lim = spe_lim(res,p_value,pcs_left)

% Control limit for SPE statistic based on J.E. Jackson and G.S. Mudholkar. 
% Control procedures for residuals associated with principal component 
% analysis. Technometrics, 21:331-349, 1979.
%
% lim = spe_lim(spe,p_value)        % complete call
%
% INPUTS:
%
% res: (NxM) Two-way residuals data matrix, N(observations) x M(variables)
%
% p_value: (1x1) p-value of the test.
%
% pcs_left: (1x1) number of components not included in the PCA model.
%
%
% OUTPUTS:
%
% lim: (1x1) control limit at a 1-p_value confidence level.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/09
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
if ndims(res)~=2, error('Incorrect number of dimensions of res.'); end;
s = size(res);
if find(s<1), error('Incorrect content of res.'); end;
if (p_value<0||p_value>1), error('Incorrect value of p_value.'); end;
if nargin < 3, pcs_left = rank(res)-1; end;
if (pcs_left<0), error('Incorrect value of pcs_left.'); end;

% Computation

lambda = sort(eig(1/(s(1)-1)*res'*res),'descend');

theta1 = sum(lambda(1:pcs_left));
theta2 = sum(lambda(1:pcs_left).^2);
theta3 = sum(lambda(1:pcs_left).^3);

h0 = 1-2*theta1*theta3/(3*theta2^2);

z = norminv(1-p_value); 
    
lim = theta1*(z*sqrt(2*theta2*h0^2)/theta1 + 1 + theta2*h0*(h0-1)/(theta1^2))^(1/h0);

h1    = z*sqrt(2*theta2*h0.^2)/theta1;
h2    = theta2*h0*(h0-1)/(theta1.^2);
rescl = theta1*(1+h1+h2).^(1/h0);

