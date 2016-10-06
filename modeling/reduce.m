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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 31/Oct/08
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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