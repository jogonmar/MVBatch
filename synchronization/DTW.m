function [sB,warp_f_ref,warping,d]= DTW(T,R,W,relaxedE,constrainedVars,band)

%  DTW function carries out the synchronization of multivariate batch trajectories using Dynamic Time Warping.
%
% CALLS:
%       [sB,warp_f_ref,warping,d]= DTW(T,R);                                 % Same weight for all variables and no global constraint defined  
%       [sB,warp_f_ref,warping,d]= DTW(T,R,W,);                              % No global constraint defined  
%       [sB,warp_f_ref,warping,d]= DTW(T,R,W,relaxedE);                      % No global constrained and relaxed end point enabled or diabled
%       [sB,warp_f_ref,warping,d]= DTW(T,R,W,false,constrainedVars);         % No global constrained and relaxed end point enabled or diabled with variables constrained
%       [sB,warp_f_ref,warping,d]= DTW(T,R,W,false,constrainedVars,band);   % Complete call
%
%
% INPUTS:
%
% T: (Kt x J) sample batch to be aligned with the reference batch. Kt is 
%     the number of measurements collected for each one of the J process
%     variables.
%
% R: (Kr x J) reference batch used to carry out the alignment. Kr
%     is the number of measurements collected for each one of the J process
%     variables.
%
% W: (JxJ) matrix containing weights to give more importance to those
%    variables are more important based on a certain criterium to achieve
%    the optimal synchronization. 
%    
%
% relaxedE: (1x1) logical value indicating whether the end point constraint is
%           relaxed or not.
%
% constrainedVar: (J x 1) weigths to give more importance to certain 
%                  process variables in batch synchronization. Note that those
%                   variables discarded must have a weight equal to 0 in the weight matrix W.
%
% band: (Kr x 2) matrix containing the upper and lower limits that define
%       the research space to estimate the local and cumulative distances, and
%       the warping path. Note that the number of points of the bands must be, at least, the same
%       as the number of sampling points of the reference batch R.
%        
% OUPUTS:
%       
% sB:  (KrxJ) synchronized batch trajectories.
%
% warp_f_ref: (Kref x I) matrix containing the warping information derived from
% batch synchronization expressed as a function of the reference batch.
%
% warping: (Kn x 2) matrix containing the warping information from the
%           offline synchronization.
%
% d: (KrxKt) local weighted distance matrix.
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
% last modification:
% October 2013: Warping information is expressed as a function of the
%               reference batch. The resulting warping profiles are equal in length
%               across batches.
% August 2013:  constrainedVars vector added as a new input in the function. It provides
%               flexiblity to constrain certain process variables in the batch
%               synchronization.
% Februrary 2013:Relaxed end point constraint is added into the algorithm.
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

 
 %% Parameters checking
 
if nargin < 2, error('Incorrect number of input paramters. Please, check the help for further details.'); end
sT = size(T);
sR = size(R);
if sT(2)~=sR(2), error('The number of process variables are not the same in T and R matrices.'); end
if iscell(R) || iscell(T), error('Reference, test or both batches are cell array. Both must be 2D matrices'); end
if ndims(sT)~=2, error('Incorrect number of dimensions of S'); end;
if ndims(sR)~=2, error('Incorrect number of dimensions of R'); end;
if nargin < 3, W = eye(sT(2),sT(2)); end
if sR(2)~= sT(2), error('The number of variables must be equal between the reference and test batch'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be squared.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
if sR(2)~=sw(1), error('The number of the weights do not match with the number of process variables.'); end
if nargin < 4 || isempty(relaxedE), relaxedE = false;end
if nargin < 5, constrainedVars = zeros(sR(2),1); end
if isempty(constrainedVars), constrainedVars = zeros(sR(2),1); end
VarNonConstr = find(constrainedVars==0);
wVars = find(diag(W)==0);
if nargin < 6, 
    if numel(VarNonConstr)~=sR(2)-numel(wVars), error('The number of constrained variables does not coincide with the null weights in the matrix W'); end
    if find(constrainedVars==0)~=find(diag(W)~=0), error('One or some variables marked as constrained does/do not have a zero weight in W'); end
    band(:,1) = ones(sT(1),1); band(:,2) = ones(sT(1), 1)* sR(1);
end
if size(band,2) ~= 2 || size(band,1) < sT(1), error('Band is not accurate for the synchronization of these multivariate trajectories'); end

% Initializing the warping vector
warp_f_ref = ones(sR(1),1).*NaN;

% Reestimation of the W matrix taking into account the constraints
W = W(VarNonConstr,VarNonConstr);


%% Initializing local distance matrix
 d = zeros(sR(1),sT(1));
 
%% DTW algorithm
 
% Iterative procedure to estimate the Mahalanobis distance between the test
% and the reference batch

 for j=1:sT(1)
    B = repmat(T(j,VarNonConstr),sR(1),1);
    d(:,j)= (diag((B - R(:,VarNonConstr))*W*(B - R(:,VarNonConstr))'))';
 end

% Cumulative distance estimation using dynamic programming based on the Sakoe-Chiba constraints.

D=ones(sR(1),sT(1))*NaN;

D(1,1)=d(1,1);

for n=2:numel(find(band(:,1)==1))
    D(1,n)=d(1,n)+D(1,n-1);
end
for m=2:band(1,2)
    D(m,1)=d(m,1)+D(m-1,1);
end

for n=2:sT(1)
    for m=max(2,band(n,1)):band(n,2)
        D(m,n)=d(m,n)+min([D(m,n-1),D(m-1,n-1),D(m-1,n)]); % Using predecessors
    end
end

n=sT(1);
m=sR(1);

endpoint = sR(1);
if relaxedE,[~,m] = min(D(:,sT(1))); endpoint = m;end

warping = [m n];

%% Step to find the optimal path with backtracking
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else
         [~,number]=min([D(m-1,n-1),D(m,n-1),D(m-1,n)]);
      switch number
      case 1
        n=n-1;
        m=m-1;
      case 2
        n=n-1;
      case 3
        m=m-1;
      end

    end
    % Storing the indeces to reconstruct the optimal path at the end of the iterative procedure
    warping = cat(1,[m n],warping);
end    


for i=1:endpoint
    warp_f_ref(i,1) = warping(find(warping(:,1)==i,1,'last'),2);
end

%% Synchronization step once the warping function has been computed

 sB = zeros(m,sR(2));
k = 1;
temp          = R(warping(1,1),:);
sB(k,:) = T(warping(1,2),:);
for i=2:length(warping(:,1))
    if warping(i,1) ~= warping(i-1,1)
       k = k+1;
       clear temp;
       temp = T(warping(i,2),:);
       sB(k,:) = T(warping(i,2),:);
    end
    if warping(i,1) == warping(i-1,1)
       temp = [temp;T(warping(i,2),:)];
       sB(k,:) = nanmean(temp);
   end
end
