function Bref = refBatch(X,method)

% refBatch function allows us to find a reference batch based on a
% criterium.
%
% INPUTS:
%
% X: (1xI) cell array containing the measurements collected for J variables at 
%    Ki different sampling times for each one of the I batches. 
%           
% method:   criterium to be followed to estimate the batch reference.
%           
% OUTPOUT:
%
% Bref:     number of the batch selected as reference following the method
%           provided.
%
% CALLS: 
%           refBatch(X,method)      % complete call
%
% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modification: 15/May/11.

%% Parameters cheching

if nargin < 2, errordlg('The number of arguments in not correct.'); end
if ~iscell(X), error('The data structure must be a cell array to contain the unsynchronized batch trajectories.'); end
    
I =length(X);
l = zeros(1,I);

for i=1:I
    l(i)=size(X{i},1);
end

switch method
    case 'median'
    [num Bref] = min(abs(l - median(l)));
    case 'average'
    [num Bref] = min(abs(l - mean(l)));
    otherwise
        errordlg('The method introduced is not correct.')
end
        