function [phases,clu]=reduce(tree)

% Computation of the final phases from the leaves of the information tree in 
% a MP model.
%
% phases=reduce(tree) % complete call
%
%
% INPUTS:
%
% tree: (n_nodesx5) complete information of the analysis with
%       the MP algorithm. Each row contains the information of a node,
%       namely [PRESS, PCs, lags, initial time, end time].
%
%
% OUTPUTS:
%
% phases: (n_phasesx5) phases of the MP model. Each row contains 
%       the information of a phase, namely [PRESS, PCs, lags, initial time,
%       end time].
%
% clu: (1xK) vector with the assignment of the K sampling times to the phases, 
%       numbered from 1 onwards.
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 31/Oct/08.

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if ndims(tree)~=2, error('Incorrect number of dimensions of tree.'); end;
s = size(tree);
if find(s(2)~=5), error('Incorrect content of tree.'); end;


% Initialization

ind=tree(:,4);
i=1;
o=1;
phases=[];
clu=[];


% Phases extraction

while i<=length(ind),
    ind2=find(ind==ind(i));
    phases=[phases;tree(ind2(end),:)];
    clu=[clu;o*ones(tree(ind2(end),5)-tree(ind2(end),4)+1,1)];
    i = ind2(length(ind2))+1;
    o = o+1;
end

if ~isinf(tree(1,3)) clu(1:tree(1,3)) = 0; end