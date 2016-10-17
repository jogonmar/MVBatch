function [matrix, rng] = scale_(rawMatrix,rng)

% scale_ function scales a set of batches to account for the differences in the units
% used to record the trajectories of the J variables. This function divides
% each variable across all batches by its average range. Such range can be
% provided or estimated.
%
% INPUT:
%
% rawMatrix: set of I batches containing the measurements collected of J
%              variables at each kn sampling time.
%
% rng:       (1xJ) vector containing the average range of the J variables
%              
% OUTPUT:
%
% matrix:    data structure containing a set of scaled batches or an unique scaled batch.
%
% rng:       (1xJ) vector containing the average range of the J variables
%
% CALLS:
%            [matrix, rng] = scale_(rawMatrix)         % call to estimate the average range-based scaled data matrix and the scaling vector.
%            [matrix, rng] = scale_(rawMatrix, rng)    % call to estimate the average range-based scaled data matrix.
%
% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modification: 15/May/11.

 %% Checking input parameters
 if iscell(rawMatrix)
     nVariables = size(rawMatrix{1},2);
 else
     nVariables = size(rawMatrix,2);
 end
  if nargin == 2 && size(rng,2)~= nVariables, errordlg('The number of weights in the scaling vector does not match with the number of the process variables'); end
  
 %% Scaling the batch process data 
if iscell(rawMatrix)

    range = zeros(size(rawMatrix,2),nVariables);
    if nargin < 2
        for i=1:size(rawMatrix,2)
      
            for j=1:nVariables

                 range(i,j) = nanmax(rawMatrix{i}(:,j)) - nanmin(rawMatrix{i}(:,j));
                 if range(i,j) == 0, range(i,j) = 0.000001;end
            end
        end
        rng = nanmean(range);
    end

    matrix = rawMatrix;


    for i=1:size(rawMatrix,2)
    for j=1:nVariables
    matrix{i}(:,j) = matrix{i}(:,j)./rng(j);
    end
    end

else
    matrix = rawMatrix;
    for j=1:size(rawMatrix,2)
    matrix(:,j) = matrix(:,j)./rng(j);
    end
end
    
