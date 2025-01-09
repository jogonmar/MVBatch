function F = create_design(levels,reps)

% Creates a balanced design matrix according to an experimental design.
%
% F = create_design(levels)   % minimum call
% F = create_design(levels,reps)   % complete call
%
%
% INPUTS:
%
% levels: {F} cell with the levels of the factors, specified as vectors.
%
% reps: [1x1] number of replicates per combination of levels in the
% factors.
%
%
% OUTPUTS:
%
% F: [NxF] design matrix.
%
%
% EXAMPLE OF USE: Two factors, with 4 and 3 levels, and 4 replicates:
%
% reps = 4;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 20/Feb/18.
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 2 || isempty(reps), reps=1; end;

%% Main code

F = ones(reps,1);
for i = length(levels):-1:1,
    Fn = repmat(F,length(levels{i}),1);
    for j = 1:length(levels{i}),
        Fn(((j-1)*size(F,1)+1):(j*size(F,1)),i) = levels{i}(j);
    end
    F = Fn;
end