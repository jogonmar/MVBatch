function y = align_IV(x,var,steps,method)

% Alignment of a batch according to an indicator variable.
%
% y = align_IV(x,var,steps)   % nearest neighbor interpolation
% y = align_IV(x,var,steps,method)   % complete call
%
%
% INPUTS:
%
% x: (KxJ) data corresponding to the J variables collected in the batch. 
%
% var: (1x1) index of the indicator variable. 
%
% steps: (1x1) number of sampling times in the aligned data.
%
% method: see methods in function interp1:
%       'nearest': nearest neighbor interpolation (default)
%       'linear': linear interpolation
%       'spline': piecewise cubic spline interpolation (SPLINE)
%       'pchip': shape-preserving piecewise cubic interpolation
%       'cubic': same as 'pchip'
%       'v5cubic': the cubic interpolation from MATLAB 5, which does not
%                   extrapolate and uses 'spline' if X is not equally
%                   spaced.
%
% OUTPUTS:
%
% y: (stepsxJ) data corresponding to the variables collected in the batch,
%   except for the indicator variable. The column 'var' of the indicator
%   variable is filled with the original sampling timing aligned according
%   to the indicator variable.
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 20/Aug/09.

% Parameters checking

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, method = 'linear'; end; 
s = size(x);
if (var<1||var>s(2)), error('Incorrect value of var.'); end;


% Delete nans
ind = find(~isnan(x(:,var)));
x2 = interp1(ind,x(ind,var),(1:s(1))'); % for interlaced nans
ind_tot = ind(1):ind(end); % for initial or final nans
if x2(1)>x2(end), x2 = -x2; end;
M = max(x2);
m = min(x2);

% Delete duplicated values
[x3,inds] = sort(x2(ind_tot));
res = [1;x3(2:end)-x3(1:end-1)];
ind_tot = ind_tot(sort(inds(find(res))));

y = interp1((x2(ind_tot)-m)/(M-m),[x(ind_tot,1:var-1) M*ind_tot'/ind_tot(end)-m x(ind_tot,var+1:end)],((0:(length(ind_tot)-1)/(steps-1):length(ind_tot)-1)/(length(ind_tot)-1))',method);

% xi = interp1(x2, ind_tot, range(1):((range(2)-range(1))/(steps-1)):range(2));
% 
% y = interp1(ind_tot, x,xi);
% y(:,6) = xi;%range(1):((range(2)-range(1))/(steps-1)):range(2);
