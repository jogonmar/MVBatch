function y = fold(x,batches,lag,method)

% Folding of batch process data.  
%
% y = fold(x,batches,lag)           % reconstruction by first value
% y = fold(x,batches,lag,method)    % complete call
%
% INPUTS:
%
% x: (NxM) two-way data matrix, N(samples) x M(variables).
%
% batches: (1x1) number of batches of the original data.
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to a sample.
%
% method: reconstruction procedure:
%   'mean': mean of all the values of a variable.
%   'first': first value (by default).
%
%
% OUTPUTS:
%
% y: (KxJxI) folded matrix, K(sampling times) x J(variables)
%   x I(batches).
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 23/Oct/06.

% Parameters checking

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, method = 'first'; end;

if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if mod(s(1),batches)~=0, error('Incorrect value of batches.'); end;
if mod(s(2),lag+1)~=0, error('Incorrect value of lag.'); end;
if lag<0, error('Incorrect value of lag.'); end;


% Folding

d(1) = s(1)/batches + lag;
d(2) = s(2)/(lag+1);
d(3) = batches;

y = ones(d(1),d(2),d(3));

for i=1:d(2),
    y(1:lag,i,:)=x(1:batches,i:d(2):(s(2)-d(2)))';
end


switch lower(method)
    
    case 'mean' % Mean
        for i=lag+1:d(1),
            yint = x((i-lag-1)*batches+1:(i-lag)*batches,lag*d(2)+(1:d(2)))'/min(lag+1,d(1)-i+1);
            for o=1:min(lag,d(1)-i),
                yint = yint + x((i-lag+o-1)*batches+1:(i-lag+o)*batches,(lag-o)*d(2)+(1:d(2)))'/min(lag+1,d(1)-i+1);
            end
            y(i,:,:) = yint;
        end
              
    case 'first' % First appearing
        for i=lag+1:d(1),
            y(i,:,:) = x((i-lag-1)*batches+1:(i-lag)*batches,lag*d(2)+(1:d(2)))';
        end
        
    otherwise
        error('Incorrect method.');
end