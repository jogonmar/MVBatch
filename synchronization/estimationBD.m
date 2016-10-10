function band  = estimationBD(warping, Kref, maxSPath)

% estimationBD estimates the bands from the warping paths obtained in off-line synchronization.
%
% INPUTS:
%
% warping: (1xI) cell array containing the warping information from the
%           off-line synchronization of the I historical batches.
%
% Kref: integer indicating the reference batch duration.
% 
% maxSPath: integer indicating the maximum batch duration of the
%       historical data set.
%
% OUTPUTS:
%
% band: (max(Ki)x 2) matrix containing the upper and lower limits which define
%       the research space to estimate the local and cumulative distances, and
%       the warping path. Note that max(Ki) is the maximum batch duration of
%       the I historical batches.
%
% 
% CALLS:
%
%       band  = estimationBD(warping, Kref, maxSPath)      % complete call
%
% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modifications 
% 17/Oct/13: The warping information of each batch is expressed as a
% function of the reference batch, yielding warping profiles of equal
% length across batches.

%% Arguments checking
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if iscell(warping) && nargin<2, error('Unexpected number of input parameters for warping profiles stored in cell arrays'); end

if iscell(warping)
    band = zeros(maxSPath,2);
    for i=1:maxSPath
        boundary = [];
        for j=1:length(warping)
            if( i<= max(warping{j}(:,2)))
            index  = find(warping{j}(:,2)==i);
            boundary = [boundary min(warping{j}(index,1)) max(warping{j}(index,1))];
            else
                boundary = [boundary Kref];
            end
        end
        band(i,:) = [min(boundary) max(boundary)];
    end
else
    band = [];
    band = [band min(warping')'];
    band = [band max(warping')'];
end




