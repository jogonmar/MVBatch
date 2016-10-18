function [sBn,warpEref,warp,band] = onSyn(Bn,Bref,band,W,zeta,rng)

% Real-time synchronization using the RGTW algorithm. 
%
% CALLS:
%
%       [sBn,warpEref,warp,band] = onSyn(Bn,Bref,band,W)               % minimum call
%       [sBn,warpEref,warp,band] = onSyn(Bn,Bref,band,W,zeta,rng)      % complete call
%
% INPUTS:
%
% Bn: (Kt x J) sample batch to be aligned with the reference batch. Kt is 
%     the number of measurements collected for each one of the J process
%     variables.
%
% Bref: (Kr x J) reference batch used to carry out the alignment. Kr
%     is the number of measurements collected for each one of the J process
%     variables.
%
% band: (Kr x 2) matrix containing the upper and lower limits that define
%       the research space to estimate the local and cumulative distances, and
%       the warping path. Note that the number of points of the bands must be, at least, the same
%       as the number of sampling points of the reference batch R.
%
% W: (JxJ) matrix containing weights to give more importance to those
%    variables are more important based on a certain criterium to achieve
%    the optimal synchronization. 
%
% zeta: windows width to use for real-time synchronization 
%      
% rng: (1xJ) vector containing the mean range of each of J variable trajectories.
%        
% OUTPUTS:
%
% sBn:  (KrefxJ) matrix containing the synchronized batch trajectories.
%
% warpEref: (Kref x 1) vector containing the warping information derived from
% batch synchronization expressed as a function of the reference batch.
%
% warp: (1 x Kref) cell array containing the warping information obtained from batch
% synchronization.
%
% band: (Kref x 2) matrix containing the updated upper and lower limits from
%        the realtime synchronization. 
%
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
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

sizeR = size(Bref); sizeS= size(Bn);

if nargin < 4, error('Incorrect number of input paramters. Please, check the help for further details.'); end
if iscell(Bref) || iscell(Bn), error('Reference, test or both batches are cell array. Both must be two-way arrays'); end
if ndims(Bn)~=2, error('Incorrect number of dimensions of S'); end;
if ndims(Bref)~=2, error('Incorrect number of dimensions of R'); end;
if sizeR(2)~= sizeS(2), error('The number of variables must be equal in the reference and test batches'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be square.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
if sizeR(2)~=sw(1), error('The number of weights do not match with the number of process variables.'); end
if zeta < -1, error('The number of window widths under study must be greater than 0'); end
if size(rng,2)~=sizeR(2), error('The number of scaled values do not match with the number of process variables.'); end
if size(band,2) ~= 2 , error('Dimension Error. Please, check the help for further details on the Band parameter.'); end

% Initialization
alignment = struct( 'nsamples',1,...
                    'temp',[],...
                    'k',1,...
                    'kprev',NaN,...
                    'aBatch',[],...
                    'kwarp',0,...
                    'warp',[],...
                    'iConst',1,...
                    'jConst',1,...
                    'consecPoi',0,...
                    'z',zeta);    

alignment.warp  = [alignment.warp; [1 1]];

Bns = scale_(Bn,rng); 

[alignment,band,d,D] = RGTW (Bns(1,:),Bref,W,band,alignment,[],[]);

for i=2:size(Bn,1)
    [alignment,band,d,D,warping]  = RGTW (Bns(1:i,:),Bref,W,band,alignment,d,D);
end

warp = alignment.warp;

% Inserting the last elements from the warping information vector

c = find(warping(:,2) == alignment.warp(size(alignment.warp,1),2),1,'last');

warp = [warp; warping(c+1:end,:)];

 a = (warp(end,1)+1:size(Bref,1))';
 b = ones(size(Bref,1)-warp(end,1),1).*warp(end,2);

warp= [warp; [a b]];

% Expressing the warping information as a function of the reference batch
warpEref = zeros(size(Bref,1),1);

for i=1:size(Bref,1)
    warpEref(i,1) = warp(find(warp(:,1)==i,1,'last'),2);
end

% Obtaining the synchronized trajectories without scaling
sBn = alignment.aBatch.*repmat(rng,size(alignment.aBatch,1),1);
                    

