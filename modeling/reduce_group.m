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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 6/Sep/07
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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