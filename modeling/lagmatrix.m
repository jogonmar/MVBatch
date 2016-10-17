function y = lagmatrix(x,lags)

% Addition of lags to two-way array
%
% y = lagmatrix(x,lags)     % complete call
%
%
% INPUTS:
%
% x: (KxJ) two-way batch data matrix, K(sampling times) x J(variables)
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
% one in the row of the unfolded matrix.
%
%
% OUTPUTS:
%
% y: ((K-lag)x(J(1+lag))) lagged matrix. [x(k) ... x(k-lag+1) x(k -lag)]
%
%
% coded by: José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2016  Technical University of Valencia, Valencia
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
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end
s = size(x);
if lags<0, error('Incorrect value of lags.'); end;
if lags>s(1)-1, lags=s(1)-1; end;

y=[];
for i=1:lags+1   
   y = [y x(i:end-lags+i-1,:)];
end