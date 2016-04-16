function mp_group2=reduce_group(mp_group)

% Erase repeated MP models in a group of models. 
%
% mp_group2=reduce_group(mp_group) % complete call
%
%
% INPUTS:
%
% mp_group: (cell) containing a set of MPPCA models.
%
%
% OUTPUTS:
%
% mp_group2: (cell) output set of MPPCA models.
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 6/Sep/07.

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;

% Main code

mp_group2={};
for i=1:length(mp_group), 
    already=false;
    for j=1:length(mp_group2), % Only add if it was not yet
        if isequal_phase(mp_group2{j}.phases,mp_group{i}.phases),
            already=true;
        end
    end
    if ~already,
        mp_group2{end+1}=mp_group{i};
    end
end