function [W,X,warping,rng,warpingOri,diffW,flag] = DTW_Kassidas(cal,ref,Wconstr,maxIter)

% Function to synchronize a set of batch trajectories against a stated
% batch (ref batch) through an iterative procedure. The synchonization
% is carried out by using the classical DTW, giving more weight to
% those variables that are more consistent batch to batch. 
% The original paper is: Kassidas et al. (1998). Synchronization of batch 
% trajectories using dynamic time warping, AIChE Journal 44:864-875.
%
% CALLS:
%        [W,X,warping,rng] = DTW_Kassidas(cal, ref)                   % minimum call
%        [W,X,warping,rng] = DTW_Kassidas(cal, ref, Wconstr,maxIter)  % complete call
%
%
% INPUTS:
%
% cal:      (1xI) cell array containing the measurements collected for J variables at 
%            Ki different sampling times for all I batches.
%       
% ref:       (KrefxI) reference batch.
%
% Wconstr: (Jx1) boolean array indicating if a specific variable is
%            considered in the synchronization (0) or not (1).
%        
% maxIter: (1x1) maximum number of iteration allowed for the Kassidas et
% al.'s synchronization. By default, maxIter is equal to 20.
%
% OUTPUTS:
%
% W:       (JxJ) weights to give more importance to those variables are more 
%          consistent batch to batch.
%
% X:       (KxJxI) data matrix containing J batch trajectories expressed at K time intervals
%          for all I batches.
%
% warping: (KrefxI) warping information expressed as a function of the
%           reference batch.
%
% rng:      (1xJ) mean range of each J variable trajectory.
%
% warpingOri:(1xI) warping information derived from the offline synchronization 
%             of I historical batches.
%
% diffW: (FxJ) array of dissimisilarities among weights across F iterations
%
%
% flag: flag to indicate if the synchronization was unexpectedly
% interrupted.
%
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
% last modification: July/11.
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

%% Parameters cheching

if nargin < 2, error('Incorrect number of input paramters. Please, check the help for further details.'), end
if ~iscell(cal), error('The first input parameter must be a cell array containing unsynchronized trajectories.'); end
nBatches = size(cal,2); nVariables = size(cal{1,1},2);
if nargin < 3 || isempty(Wconstr) , Wconstr = zeros(size(cal{1,1},2),1); end
if find(Wconstr==0) < 1, error('If all the process variables are constrained, no synchronization can be performed.');end
if nargin < 4, maxIter = 20; end
if isempty(maxIter), maxIter = 20; end

%% Initialization
calsy = cell(1, nBatches);
warping = zeros(size(ref,1), nBatches);
warpingOri = cell(1, nBatches);
X = zeros(size(ref,1),nVariables,nBatches);
W = zeros(nVariables,1);
diffW = [];
Bmean = [];
flag = true;
iter = 0;
% Updtate of the weighting matrix based on the variable constraints
VarConstr = find(Wconstr==0);
W(VarConstr) = nVariables/numel(VarConstr);
W = diag(W);


%% Step A: Scaling
[calsc,rng] = scale_(cal); 
refs = scale_(ref,rng);

%% Step B: Synchronization

 h = waitbar(0/nBatches,sprintf('Synchronizing batches - Iteration %d',iter+1),'Name','Synchronization by Kassidas et al. DTW algorithm');

 
while(flag || iter <=4)
    
    if iter < 3
        Bref = refs;
    else
        Bref = Bmean;
    end
    % Step 1: Apply the DTW/synchronization method between Bi i=1,...,I y Bref
    for i=1:nBatches          
         % Check for Cancel button press
        if getappdata(h,'canceling')
            delete(h);
            return;
        end
        waitbar(i/nBatches,h,sprintf('Synchronizing batch #%d - Iteration %d',i,iter+1),'Name','Synchronization by Kassidas et al. DTW algorithm');
        [calsy{i},warping(:,i),warpingOri{i}] = DTW(calsc{i},Bref,W,0,Wconstr);
    end

    % Step 2: Compute the average trajectory 
    Add = zeros(size(refs,1),nVariables);
    for i=1:nBatches
        Add = Add + calsy{i}(:,:);
    end
    Bmean = Add ./ nBatches;

    % Step 3: For each variable, compute the sum of squared desviations from Bmean
    Wnew = diag(zeros(nVariables,1));
    for j=1:numel(VarConstr)
        weight = 0;
        for i=1:nBatches
            weight = weight + (nansum((calsy{i}(:,VarConstr(j)) - Bmean(:,VarConstr(j))).^2));
        end     
        Wnew(VarConstr(j),VarConstr(j)) = 1/weight;      
        if weight == 0, Wnew(VarConstr(j),VarConstr(j)) = 0.0001;end
    end
    % Step 4: Normalize W so that the sum of the weights is equal to the number of
    %         variables.

    Wnew = diag((diag(Wnew)./sum(diag(Wnew)))*nVariables);    
    diffW = [diffW; diag(abs(Wnew) - abs(W))'];
    if  isempty(find(diffW(end,:) > 0.01)) || iter == maxIter
        flag = false;
    end 
    W = Wnew;
    iter = iter + 1;
    
end

for i=1:nBatches            
    [calsy{i},warping(:,i),warpingOri{i}] = DTW(calsc{i},refs,W,0,Wconstr);
end

close(h)

%% Obtaining the synchronized trajectories without scaling

for i=1:nBatches
    X(:,:,i)=calsy{i}.*repmat(rng,size(calsy{i},1),1);
end

%% Setting the remaining outputs
W = diag(W);
