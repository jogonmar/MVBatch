function [sB, warp, warping, d]= DTW(T,R,W,relaxedE,constrainedVars,band)

%  DTW function carries out the synchronization of multivariate batch trajectories using Dynamic Time Warping.
%
% INPUTS:
%
% T: (Kt x J) reference batch which is used to carry out the alignment. Kt
%     are the number of measurements collected for each one of the J process
%     variables.
%
% R: (Kr x J) sample batch which will be align against to the reference
%     one. Kr are the number of measurements collected for each one of the J process
%     variables.
%
% W: (JxJ) matrix containing weights to give more importance to those
%    variables are more proper to achieve the synchronization based on a
%    criterium.
%
% relaxedE: logical value indicating whether the end point constraint is
% relaxed or not.
%
% constrainedVar: (J x 1) vector indicating the process variables that will
% have a certain importance in the synchronization. Note that those
% variables discarded must have a weight 0 in the weight matrix W.
%
% band: (Kr x 2) matrix containing the upper and lower limits which define
%       the research space to estimate the local and cumulative distances, and
%       the warping path. Note that the number of points of the bands must be, at least, the same
%       as the number of sampling points of the reference batch R.
%        
% OUPUTS:
%       
% sB:  (KrxJ) matrix containing the synchronized batch trajectories.
%
% warp: (Kref x I) matrix containing the warping information derived from
% batch synchronization expressed as a function of the reference.
%
% warping: (Kn x 2) matrix containing the warping information from the
%           off-line synchronization.
%
% dist: (KrxKt) local weighted distance matrix
%
% distm: mean cummulated distance independent of the number of synchronization steps carried out.
%
% CALLS:
%       [sB,warping,distm]= DTW(T,R);                                 % Same weight for all variables and no global constraint defined  
%       [sB,warping,distm]= DTW(T,R,W,);                              % No global constraint defined  
%       [sB,warping,distm]= DTW(T,R,W,relaxedE);                      % No global constrained and relaxed end point enabled or diabled
%       [sB,warping,distm]= DTW(T,R,W,false,constrainedVars);         % No global constrained and relaxed end point enabled or diabled with variables constrained
%       [sB,warping,distm]= DTW(T,R,W,false,constrainedVars, band);   % Complete call
%
% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modification: 
% 17/Aug/13: constrainedVars vector added as new input in the function. It provides
%            flexiblity to constrain certain process variables in the batch
%            synchronization.
% 07/Feb/13: Relaxed end point constraint option is incorporate in the algorithm.
% 01/Jan/13: An output parameter corresponding to the average cummulated distance independent of the number of
% synchronization steps has been added.
% 17/Oct/13: Warping information is expressed as a function of the
% reference batch. The resulting warping profiles are equal in length
% across batches.
 


 %% Parameters checking
 
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
warp = ones(sR(1),1).*NaN;

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
if relaxedE,[value m] = min(D(:,sT(1))); endpoint = m;end

warping = [m n];

%% Step to find the optimal path with backtracking
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else
         [values,number]=min([D(m-1,n-1),D(m,n-1),D(m-1,n)]);
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
    % Storing the indeces to rebuild the optimal path at the end of the iterative procedure
    warping = cat(1,[m n],warping);
end    


for i=1:endpoint
    warp(i,1) = warping(find(warping(:,1)==i,1,'last'),2);
end

%% Synchronization step once the warping function has been reached

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
%        ind=1;
%        if size(temp,1)>1, ind=2;end
%        sB(k,:) = temp(ind,:);
   end
end
