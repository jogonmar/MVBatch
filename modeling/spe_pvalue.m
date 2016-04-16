function p_value = spe_pvalue(res,lim,pcs_left)

% Inverse of the control limit for SPE statistic based on J.E. Jackson and 
% G.S. Mudholkar. Control procedures for residuals associated with principal 
% component analysis. Technometrics, 21:331-349, 1979.
%
% p_value = spe_pvalue(res,lim,pcs_left)       % complete call
%
% INPUTS:
%
% res: (NxM) Two-way residuals data matrix, N(observations) x M(variables)
%
% lim: (1x1) control limit at a 1-p_value confidence level.
%
% pcs_left: (1x1) number of components not included in the PCA model.
%
%
% OUTPUTS:
%
% p_value: (1x1) p-value of the test.
%
%
% codified by: José Camacho Páez.
% last modification: 05/May/09.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(res)~=2, error('Incorrect number of dimensions of res.'); end;
s = size(res);
if find(s<1), error('Incorrect content of res.'); end;
if (lim<0), error('Incorrect value of lim.'); end;
if nargin < 3, pcs_left = rank(res)-1; end;
if (pcs_left<0), error('Incorrect value of pcs_left.'); end;

% Computation

lambda = eig(1/(s(1)-1)*res'*res);

theta1 = sum(lambda(1:pcs_left));
theta2 = sum(lambda(1:pcs_left).^2);
theta3 = sum(lambda(1:pcs_left).^3);

h0 = 1-2*theta1*theta3/(3*theta2^2);

p_value = 1-normcdf(theta1*((lim/theta1)^h0 - 1 - theta2*h0*(h0-1)/(theta1^2))/sqrt(2*theta2*h0^2));
 


