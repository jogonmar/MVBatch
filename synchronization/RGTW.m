function [alignment,band,d,D,warping] = RGTW (S,R,W,band,alignment,d,D)

% RGTW algorithm to sincronize two batch trajectories in online
% applications.
%
% The original work is: 
% [1] González et al. Real-time synchronization of batch trajectories for on-line 
% multivariate statistical process control using Dynamic Time Warping, Chemometrics 
% and Intelligent Laboratory Systems, 105 (2011) 195-206).
% [2] González-Martinez, J.M. Advances on bilinear modeling of biochemical
% batch processes (2015). PhD thesis, DOI: 10.4995/Thesis/10251/55684.
%
% INPUTS:
%
% S: (KsxJ) pre-processed sample batch which be aligned with the reference batch, 
%     where Ks is the sampling time and J the number of variables.
%
% R: (KrxJ) pre-processed reference batch which is used to perform synchronization,
%    where Kr is the sampling time and J the number of variables.
%
% W: (JxJ) weight matrix used used in the local distance computation in order to give more
%    importance to certain variables based on a criterium.
%
% band: (max(Kn) x 2) upper and lower bands.
%
% alignment: data structure containing information belonging to the online synchronization 
%
% d: local distance matrix which stores distances from predecessors points
%
% D: cumulative distance matrix which stores cumulatived distance for
%    both trajectories in the sampling time k.
%
%
% OUTPUTS:
%
% alignment: data structure containing information belonging to the online synchronization
%
% band: (max(Kn) x 2) updated upper and lower bands.
%
% d: local distance matrix which stores weighted distances in each k
%    sampling point.
%
% D: acumulative distance matrix which stores cumulatived distance for
%    both trajectories in the sampling time point k.
%
% warping: array with the warping Kref from the reference and test batches.
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
% last modification: 11/Nov/09 version 2.0 (modification of the bands)
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


%% Parameters checking and settings

sizeR = size(R); sizeS= size(S);
warping=[];

if nargin < 5, error('Incorrect number of input paramters. Please, check the help for further details.'); end
if iscell(R) || iscell(S), error('Reference, test or both batches are cell array. Both must be 2D matrices'); end
if ndims(S)~=2, error('Incorrect number of dimensions of S'); end;
if ndims(R)~=2, error('Incorrect number of dimensions of R'); end;
if sizeR(2)~= sizeS(2), error('The number of variables must be equal in the reference and test batches'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be square.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
if sizeR(2)~=sw(1), error('The number of weights do not match with the number of process variables.'); end
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

% Cumulative distance estimation at the current sampling time point k using
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
 
    % If there are warping coordenates from last samples, we arrange them for further treatment in this iteration.    
    ini = size(alignment.warp,1);
    
    indxBefore =  find(warping(:,2)==alignment.nsamples-alignment.z+1);
    if numel(indxBefore) > 1
         alignment.warp = [ alignment.warp; warping(indxBefore(2:length(indxBefore)),:)];
    end
    
    % Getting new sample to be projected onto the latent subspace. We
    % arrange into the warping vector
    
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
        % This case is executed in the first iteration since it is the unique case different from others.
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
    
