function y = unfold(x,lag)

% Unfolding of batch process data.  
%
% y = unfold(x)         % for variable-wise unfolding
% y = unfold(x,Inf)     % for batch-wise unfolding
% y = unfold(x,lag)     % complete call
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
% one in the row of the unfolded matrix (0 by default).
%
%
% OUTPUTS:
%
% y: ((K-lag)x(J(1+lag))) unfolded matrix. [x(k -lag) x(k-lag+1) ... x(k)]
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 22/Oct/06.

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if nargin < 2, lag = 0; end;

if ndims(x)~=3, 
    if ndims(x)~=2, 
        error('Incorrect number of dimensions of x.');
    else
        s = [size(x),1];
        x2 = ones(s);
        x2(:,:,1) = x;
    end
else
    s = size(x);
end

if find(s<1), error('Incorrect content of x.'); end;
if lag<0, error('Incorrect value of lag.'); end;

if lag>s(1)-1, lag=s(1)-1; end;

% Unfolding

y = ones(s(3)*(s(1)-lag),s(2)*(lag+1));
yint = ones(s(2),s(3));

for ind=0:lag,
    for o=(ind+1):(s(1)-lag+ind),
        yint(:,:) = x(o,:,:);
        y((o-ind-1)*s(3)+1:(o-ind)*s(3),(ind*s(2)+1):(s(2)*(ind+1))) = yint';
    end
end
