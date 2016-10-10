function [W,X,warp,rng,warping]=DTW_Ramaker(cal,ref,Wconstr)

% Synchronization of batch trajectories giving more importance to those 
% variables containing important warping information. 
% The original paper is: 
% [1] Henk-Jan Ramaker,  Eric N.M. van Sprang, Johan A. Westerhuis, Age K. Smilde (2003).
% Dynamic Time Warping of spectroscopic batch data, Analytica Chimica Acta 498: 133-155.
%
% CALLS:
%
%        [W,X,warping,rng]=DTW_Ramaker(cal,ref)               % minimum call
%        [W,X,warping,rng]=DTW_Ramaker(cal,ref,Wconstr)       % complete call
%
% INPUTS:
%
% cal: (1xI) cell array containing the measurements collected for J variables at 
%       Ki different sampling times for each one of the I batches.
%       
% ref: (KrefxI) refernce batch.
%
% OUTPUTS:
%
% W: (JxJ) matrix containing weights to give more importance to those
%    variables are more consistent batch to batch.
%
% X:  (KxJxI) data matrix containing the J batch trajectories, whose duration is equal to
%      K, of each one the I batches.
%
% warp: (Kref x I) matrix containing the warping information derived from
% batch synchronization.
%
% rng: (1xJ) vector containing the mean range of each one the J
%      trajectories.
%
% warping: (1xI) cell array containing the warping information from the
%           off-line synchronization of the I historical batches.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)        
% last modification: 
%
% October 2013: Warping information is expressed as a function of the
% ref batch. The resulting warping profiles are equal in length
% across batches.
%
% August 2011: Synchronization is carried out for all the batches except  for the
% ref batch. Also, the output matrix contains the synchronized test batches adding 
% the range removed in the preprocessing step.
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

%% Arguments checking

if nargin < 2, error('Number of arguments are incorrect. Please check it up.'); end
if ~iscell(cal), error('The data set has to be a cell array to save possible uneven batches.'); end
if isempty(cal), error('Error in arguments. The data set does not have to be empty'); end
if nargin < 3, Wconstr = zeros(size(cal{1,1},2),1); end
if numel(find(Wconstr==0)) + numel(find(Wconstr==1)) ~= size(cal{1},2), error('The constraints must be binary numbers. Zero values stand for variables contrainted not to have weight and one values for variables that are taken into account in the batch synchronization'); end
if size(Wconstr,1) ~= size(cal{1},2), error('The number of constraints differs from the number of variables. Please, introduce as many constraints as process variables.'); end


%% Initialization
nVariables = size(cal{1},2);
nBatches = length(cal);
X = zeros(size(ref,1),nVariables,nBatches);
VarNonConstr = find(Wconstr==0);

%% Scaling
[calsc,rng] = scale_(cal); 
Bref = scale_(ref,rng);

%% Iterative procedure to estimate the Mahalanobis distance between the sample
%% and the ref batch

% Wight matrix used in the local distance computation in order to give more
% importance to some variables, in this case, those containing important
% warping information.

W = ones(nVariables,1);
W(find(Wconstr==1),1) = 0;
weights = zeros(nBatches, nVariables);
d = cell(nBatches, nVariables);

flag = true;
while (flag)  
    aBatch = zeros(size(Bref,1),nVariables,nBatches);
    warping = cell(1,nBatches);
    warp = zeros(size(Bref,1),nBatches);

    % Synchronization using the matrix W calculated by Ramaker's approach
    for i=1:nBatches
        [aBatch(:,:,i),warp(:,i), warping{i}]= DTW(calsc{i},Bref,diag(W),0,Wconstr);
    end

    % Estimation of local distance matrix for each of the J variable for all batches.

    for i=1:nBatches
        for j=1:length(VarNonConstr)
            d{i,VarNonConstr(j)} = (repmat(Bref(:,VarNonConstr(j)),1,size(calsc{i},1)) - repmat(calsc{i}(:,VarNonConstr(j))',size(Bref,1),1)).^2;
            
            N1=zeros(size(d{i,VarNonConstr(j)}));
            
            % Estimating the average local distance corresponding to the coordinates lying outside the optimal path and those lying on the optimal path.     
            N1(sub2ind(size(d{i,VarNonConstr(j)}),warping{i}(:,1),warping{i}(:,2)))=1;
            N0 = logical(not(N1));
            N1=logical(N1);
            value =mean(d{i,VarNonConstr(j)}(N0))/mean(d{i,VarNonConstr(j)}(N1));
            if value == Inf, value=NaN;end
            weights(i,VarNonConstr(j)) = value;
        end
    end

    % Run the complete data set to find the new weights
    Wnew = nanmean(weights);
    Wnew = ((Wnew(:)./sum(Wnew))*size(Wnew,2));

    if  not(isempty(find((abs(Wnew - W)./  abs(W)) > 0.002)))    
        W = Wnew;
    else
        flag = false;
    end 
end

%% Obtaining the synchronized trajectories without scaling
for i=1:nBatches
    for j=1:nVariables
        X(:,j,i)=aBatch(:,j,i).*rng(j);
    end
end
