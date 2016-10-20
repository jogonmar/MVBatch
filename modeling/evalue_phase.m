function [exists,phases2,cumpress,ind1,ind2] = evalue_phase(phases,ini,fin)

% Evalue if a certain interval is represented by a phase or phases in a set 
% of phases.
%
% [exists,phases2,cumpress,ind1,ind2] = evalue_phase(phases,ini,fin) %
%   complete call
%
% INPUTS:
%
% phases: (n_phasesx5) set of disjoints phases. Each row contains the 
%   information of a phase, namely [PRESS, PCs, lags, initial time, end time].
%
% ini: (1x1) initial sampling time of the interval.
%
% fin: (1x1) final sampling time of the interval.
%
%
% OUTPUTS:
%
% exists: (boolean) 'true' if the phases/phase exist/s, 'false' otherwise.
%
% phases2: (n2_phasesx5) phase or phases representing the interval.
%
% cumpress: (1x1) cummulative PRESS of the interval.
%
% ind1: (1x1) index of the initial set of phases where the interval starts.
%
% ind2: (1x1) index of the initial set of phases where the interval ends. 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Nov/06
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

if nargin < 3, error('Error in the number of arguments.'); end;

if ndims(phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp = size(phases);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of phases.'); end;
if find(phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(phases(:,4:5)<1), error('Incorrect content of phases.'); end;

if ini > fin , error('Error in the specification of the interval.'); end;


% Search

ind1=find(phases(:,4)==ini);
ind2=find(phases(:,5)==fin);

l1=length(ind1);
l2=length(ind2);
    
if l1>1 || l2>1,
    error('Error in the content of phase.');
else
    if l1>0 && l2>0,

        exists=true;
        phases2=phases(ind1(l1):ind2(l2),:);    
        cumpress = sum(phases2(:,1));

    else
        exists=false;
        phases2=[];
        cumpress = 0;
    end
end