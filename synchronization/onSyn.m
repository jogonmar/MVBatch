function [sBn,warpEref,warp, band] = onSyn(Bn,Bref,band,W,zeta,rng)

% onSyn function performs the real-time synchronization using the RGTW
% algorithm. 

% INPUTS:
%
% Bn: (Kt x J) reference batch which will be align against to the reference
%     one. Kr are the number of measurements collected for each one of the J process
%     variables.
%
% Bref: (Kref x J) reference batch which is used to carry out the alignment. Kt
%     are the number of measurements collected for each one of the J process
%     variables.
%
% band: (Kref x 2) matrix containing the upper and lower limits which define
%       the research space to estimate the local and cumulative distances, and
%       the warping path. Note that the number of points of the bands must be, at least, the same
%       as the number of sampling points of the reference batch R.
%
% W: (JxJ) matrix containing weights to give more importance to those
%    variables are more proper to achieve the synchronization based on a
%    criterium.
%
% zeta: windows width to be used for real-time synchronization 
%      
% rng: (1xJ) vector containing the mean range of each one the J
%      trajectories.
%        
% OUTPUTS:
%
% sBn:  (KrefxJ) matrix containing the synchronized batch trajectories.
%
%
% warpEref: (Kref x 1) vector containing the warping information derived from
% batch synchronization expressed as a function of the reference batch.
%
% warp: (1 x Kref) cell array containing the warpinf information from batch
% synchronization.
%
% band: (Kref x 2) matrix containing the updated upper and lower limits from
%        the real-time synchronization. 
%
% 
% CALLS:
%
%       [sBn warp band] = onSyn(Bn,Bref, band,W, zeta, rng)      % complete call
%
% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modification: 15/May/11.

%% Parameters checking

sizeR = size(Bref); sizeS= size(Bn);

if nargin < 4, error('Incorrect number of arguments.'); end
if iscell(Bref) || iscell(Bn), error('Reference, test or both batches are cell array. Both must be 2D matrices'); end
if ndims(Bn)~=2, error('Incorrect number of dimensions of S'); end;
if ndims(Bref)~=2, error('Incorrect number of dimensions of R'); end;
if sizeR(2)~= sizeS(2), error('The number of variables must be equal between the reference and test batch'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be a squared one.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
if sizeR(2)~=sw(1), error('The number of the weights do not match with the number of process variables.'); end
if zeta < -1, error('The number of window width to be studied must be greater than 0'); end
if size(rng,2)~=sizeR(2), error('The number of scaled values do not match with the number of process variables.'); end
if size(band,2) ~= 2 , error('Band is not accurate for the synchronization of these multivariate trajectories'); end


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
                    

