function [Xs] = stageBreaker(X,stages)

% stageBreaker split up the three-away into S three-away array contaning
% the process data for each one of the stages. Note that 
%
% INPUTS:
%
% X: (Ix1) cell array containing the measurements of J variables measured at 
%       Ks different sampling time points in all I batches.
%            
% stages: (S x 1) array contaning the stage ID for the processnumber of stages the process has
%
% OUTPUTS:
%
% X: (Ix1) cell array containing I cell arrays (as many as stages) with the measurements of J variables measured at 
%       Ks different sampling time points for all the I batches stages.
%
% codified by: J.M. Gonzalez-Martinez.
% version: 1.0
% last modification: May/11.

%% Parameters cheching

if nargin < 1, error('The number of argument is not correct.'); end
if ~iscell(X), error('The X matrix must be a cell array containing batch trajectories with non-constant lenghts'); end

nBatches = length(X);
nStages = length(stages);
Xs = cell(nStages,1);


for s=1:nStages
    Xs{s} = cell(nBatches,1);
    for i=1:nBatches
        Xs{s}{i} = X{i}(find(ceil(X{i}(:,1))==stages(s)),2:end);
    end
end