function [matrix,rng] = scale_(calibration,rng)

% Scale a set of batches to account for differences in magnitude across
% variables.
%
% INPUT:
%
% calibration: set of I batches containing the measurements collected of J
%              variables at each kn sampling time. To compute the ranges and
%              scale batch data, the sata set must be provided in a cell 
%              array. If a two-way array needs to be scaled with given ranges,
%              the sata set must be provided in a (KxJ) matrix.
%
% rng:       (1xJ) vector containing the average range of J variables
%              
% OUTPUT:
%
% matrix:    data structure containing a set of scaled batches or an unique scaled batch.
%
% rng:       (1xJ) vector containing the average range of the J variables
%
% CALLS:
%            [matrix, rng] = scale_(calibration)         % call to estimate the variables ranges and scale data matrix with these values
%            [matrix, rng] = scale_(calibration, rng)    % call to scale a data matrix using the variable ranges
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
% last modification: May/11.
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

%% Parameters checking and initialization
 if iscell(calibration)
     nVariables = size(calibration{1},2);
 else
     nVariables = size(calibration,2);
 end
  if nargin == 2 && size(rng,2)~= nVariables, errordlg('The number of weights in the scaling vector does not match with the number of process variables'); end
  
 %% Scaling batch process data 
if iscell(calibration)

    range = zeros(size(calibration,2),nVariables);
    if nargin < 2
        for i=1:size(calibration,2)
      
            for j=1:nVariables

                 range(i,j) = nanmax(calibration{i}(:,j)) - nanmin(calibration{i}(:,j));
                 if range(i,j) == 0, range(i,j) = 0.000001;end
            end
        end
        rng = nanmean(range);
    end

    matrix = calibration;


    for i=1:size(calibration,2)
    for j=1:nVariables
    matrix{i}(:,j) = matrix{i}(:,j)./rng(j);
    end
    end

else
    matrix = calibration;
    for j=1:size(calibration,2)
    matrix(:,j) = matrix(:,j)./rng(j);
    end
end
    
