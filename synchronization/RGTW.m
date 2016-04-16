function [alignment band d D warping] = RGTW (S,R,W,band,alignment,d,D)

% RGTW algorithm to sincronize two batch trajectories in on-line applications
%
% INPUTS:
%
% S: (KsxJ) pre-processed sample batch which be align against to a
%    reference batch, where Ks is the sampling time and J the number of variables.
%
% R: (KrxJ) pre-processed reference batch which is used to achieve the alignment,
%           where Kr is the sampling time and J the number of variables.
%
% W: weight matrix used used in the local distance computation in order to give more
%    importance to some variables based on a criterium.
%
% band: (max(Kn) x 2) matrix containing the upper and lower bands.
%
% alignment: data structure containing information belonging to the on-line synchronization 
%
% d: local distance matrix which stores distances from predecessors points
%
% D: acumulative distance matrix which stores cumulatived distance for
%    both trajectories in the sampling time k.
%
%
% OUTPUTS:
%
% alignment: data structure containing information belonging to the on-line synchronization
%
% band: (max(Kn) x 2) matrix containing the updated upper and lower bands.
%
% d: local distance matrix which stores weighted distances at each k
%    sampling point.
%
% D: acumulative distance matrix which stores cumulatived distance for
%    both trajectories in the sampling time k.
%
% warping: array with the warping Kref from the reference and test batches.
%
% Function based on the algorithm proposed by González et al. [Real-time 
% synchronization of batch trajectories for on-line multivariate statistical
% process control using Dynamic Time Warping, Chemometrics and Intelligent
% Laboratory Systems, 105 (2011) 195-206)].
%
% codified by: Jose M. Gonzalez-Martinez
% last modification: 11/Nov/09 version 2.0 (modification of the bands)

%% Parameters checking and settings

sizeR = size(R); sizeS= size(S);
warping=[];

if nargin < 5, error('Error in the number of arguments.'); end
if iscell(R) || iscell(S), error('Reference, test or both batches are cell array. Both must be 2D matrices'); end
if ndims(S)~=2, error('Incorrect number of dimensions of S'); end;
if ndims(R)~=2, error('Incorrect number of dimensions of R'); end;
if sizeR(2)~= sizeS(2), error('The number of variables must be equal between the reference and test batch'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be a squared one.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
if sizeR(2)~=sw(1), error('The number of the weights do not match with the number of process variables.'); end
if size(band,2)~= 2, error('Dimensions of the array not expected.'); end
if ~isstruct(alignment), error('The parameter "alignment" must be a struct.'); end
if nargin < 6, d = [];end
if nargin < 7, D = [];end


%% Local and cumulative local distante matrix estimation

% Mahalanobis distance estimation between the reference and test batch at
% the current k sampling point.

localdist = ones(1,sizeR(1))*NaN;

for j=band(alignment.nsamples,1):band(alignment.nsamples,2)
    localdist(j) = (S(alignment.nsamples,:) - R(j,:))*W*(S(alignment.nsamples,:) - R(j,:))';
end

% Adding the local distances from the last sampling point to the local distances matrix
    d = [d localdist'];

% Cumulative distance estimation at the current k sampling time using
% dynamic programming based on the Sakoe-Chiba constraints.

    D = [D (ones(1,sizeR(1))*NaN)']; 

if alignment.nsamples ~= 1 
    if (band(alignment.nsamples,1)==1)
    D(1,alignment.nsamples)=d(1,alignment.nsamples)+D(1,alignment.nsamples-1);
    end
    for i=max(2,band(alignment.nsamples,1)):band(alignment.nsamples,2)
        D(i,alignment.nsamples)=d(i,alignment.nsamples)+min([D(i-1,alignment.nsamples),D(i-1,alignment.nsamples-1),D(i,alignment.nsamples-1)]); % Using predecessors
    end
else 
    D(1,1) = d(1,1);
    for i=2:band(alignment.nsamples,2)
    D(i,1)=d(i,1)+D(i-1,1);
    end
end

if(alignment.nsamples >= alignment.z)

    [~, indeces]=sort(D(:,alignment.nsamples));
    z=1;

    if (size(alignment.warp,1) > 1)
    low = alignment.warp(size(alignment.warp,1),1);
    else
        low=1;
    end
    i = indeces(z);

    while(i < low) 
       z = z + 1;
       i = indeces(z);
    end

    j=alignment.nsamples;  
     
    band = updateBound(i, alignment.nsamples, band,sizeR(1));
    
%% Step to find the optimal path with backtracking
    
warping=[];
warping = [warping; i j]; 

while (i+j)~=(alignment.iConst+alignment.jConst)
    if (i-1) == alignment.iConst-1
        j=j-1;
    elseif (j-1) == alignment.jConst-1
        i=i-1;
    else
      [~,npredecessor]=min([D(i-1,j-1), D(i-1,j), D(i,j-1)]);
      switch npredecessor
      case 1
        i=i-1;
        j=j-1;    
      case 2
        i=i-1;
      case 3
        j=j-1;
      end
    end
 
    % Storing indeces in order to rebuild the optimal path at the end of the iterative procedure
    warping = [[i j]; warping];
 
end

% Moving monitoring window if the number of samples is greater than the window size. It means we save the suboptimal path (up to nsamples-wsize) 
 
    % If there are warping coordenates from last samples, we collected them to
    % be treated in this iteration. 
    
    ini = size(alignment.warp,1);
    
    indxBefore =  find(warping(:,2)==alignment.nsamples-alignment.z+1);
    if numel(indxBefore) > 1
         alignment.warp = [ alignment.warp; warping(indxBefore(2:length(indxBefore)),:)];
    end
    
    % Getting new sample to be projected onto the latent subspace. We collect it
    % in the warping vector
    
    if alignment.z > 1
        indexWithin = find(warping(:,2)==alignment.nsamples-alignment.z+2);
        alignment.warp = [ alignment.warp; warping(indexWithin,:)];
        alignment.warpRemain = warping(indexWithin+1:size(warping,1),:);
    else
        alignment.warp = [ alignment.warp; warping(2:size(warping,1),:)];
    end
    alignment.iConst = alignment.warp(size(alignment.warp,1),1);
    alignment.jConst = alignment.warp(size(alignment.warp,1),2);
    alignment.kwarp = alignment.kwarp + 1;
    fin = size(alignment.warp,1);
    alignment.kprev = alignment.k;
    
    if (alignment.kwarp == 1)      
        % This case is executed in the first iteration since it is the unique case is different from the other ones.
        alignment.temp                  = S(alignment.warp(1,2),:);
        alignment.aBatch(alignment.k,:) = S(alignment.warp(1,2),:);
    end
      
    % Asymmetric step
    for i=ini+1:fin

        if alignment.warp(i,1) ~= alignment.warp(i-1,1)
           alignment.k = alignment.k+1;
           alignment.temp = [];
           alignment.temp = S(alignment.warp(i,2),:);
           alignment.aBatch(alignment.k,:) = S(alignment.warp(i,2),:);

        elseif alignment.warp(i,1) == alignment.warp(i-1,1)
           alignment.temp = [alignment.temp;S(alignment.warp(i,2),:)];
           alignment.aBatch(alignment.k,:) = nanmean(alignment.temp);  
       end
    end
end
    alignment.nsamples = alignment.nsamples + 1;
    
