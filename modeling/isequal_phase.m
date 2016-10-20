function eq = isequal_phase(phases,phases2)

% Checks if two phases are equal. 
%
% eq = isequal_phase(phase,phase2) % complete call
%
%
% INPUTS:
%
% phases: (n_phasesx5) first set of disjoints phases modelling an interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
% 
% phasse2: (n2_phasesx5) second set of disjoints phases modelling an interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
%
%
% OUTPUTS:
%
% eq: (boolean) 'true' if phases are equal and 'false' otherwise.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 20/Nov/06
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

if nargin < 2, error('Error in the number of arguments.'); end;

if ndims(phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp1 = size(phases);
if (sp1(1)<1||sp1(2)~=5), error('Incorrect content of phases.'); end;
if find(phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(phases(:,4:5)<1), error('Incorrect content of phases.'); end;

if ndims(phases2)~=2, error('Incorrect number of dimensions of phases2.'); end;
sp2 = size(phases2);
if (sp2(1)<1||sp2(2)~=5), error('Incorrect content of phases2.'); end;
if find(phases2(:,1:3)<0), error('Incorrect content of phases2.'); end;
if find(phases2(:,4:5)<1), error('Incorrect content of phases2.'); end;


% Comparison

if sp1(1)==sp2(1),
    if ~isempty(find(isinf(phases(:,3)))) || ~isempty(find(isinf(phases2(:,3)))),
        value = sum(sum((phases(:,[2 4:5])-phases2(:,[2 4:5])).^2));
    else
        value = sum(sum((phases(:,2:5)-phases2(:,2:5)).^2));
    end
end

if sp1(1)==sp2(1) && value==0,
    eq=true;
else
    eq=false;
end