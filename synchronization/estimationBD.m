function band  = estimationBD(warping,kref,max_duration)

% Estimate the bands wrapping the warping profiles obtained from batch synchronization.
%
% CALLS:
%
% band  = estimationBD(warping)                        % complete call if
% warping is a two-way array
% band  = estimationBD(warping,kref,max_duration)      % complete call
% if warping is a cell array
%
% INPUTS:
%
% warping:  (1xI) cell array containing warping information from the
%           offline synchronization of I historical batches or (KxI)
%           warping profiles.
%
% kref:     (1x1) reference batch duration.
% 
% max_duration: (1x1) maximum batch duration of the historical data set.
%
% OUTPUTS:
%
% band:     (max(Ki)x 2) matrix containing the upper and lower limits that define
%           the research space to estimate the local and cumulative distances, and
%           the warping path. Note that max(Ki) is the maximum batch duration of
%           I historical batches. If a two-way array is provided as an
%           input for the warping profiles, this output will be a (Kx2)
%           two-way array.
%
% 
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)
% last modification: 
% 17/Oct/13: The warping information of each batch is expressed as a
% function of the reference batch, leading to warping profiles of equal
% length across batches.
%
% Copyright (C) 2016  José M. Gonzalez-Martinez
% Copyright (C) 2016  Technical University of Valencia, Valencia
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

% Parameter checking
if nargin<1, error('Incorrect number of input paramters. Please, check the help for further details.'); end
if nargin<2 && iscell(warping), error('Unexpected parameters. Only providing the warping profiles, the single input parameter must be a two-way array.');end
if iscell(warping) && nargin <3, error('Incorrect number of input paramters. Please, check the help for further details.'); end
if iscell(warping) && kref<1, error('Incorrect batch duration.'); end
if iscell(warping) && max_duration < 1, error('Incorrect batch duration.'); end

% Computation
if iscell(warping)
    band = zeros(max_duration,2);
    for i=1:max_duration
        boundary = [];
        for j=1:length(warping)
            if( i<= max(warping{j}(:,2)))
            index  = find(warping{j}(:,2)==i);
            boundary = [boundary min(warping{j}(index,1)) max(warping{j}(index,1))];
            else
                boundary = [boundary kref];
            end
        end
        band(i,:) = [min(boundary) max(boundary)];
    end
else
    band = [];
    band = [band min(warping')'];
    band = [band max(warping')'];
end




