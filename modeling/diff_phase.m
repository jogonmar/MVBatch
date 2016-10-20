function [diff_phas,diff_phas2] = isequal_phase(phases,phases2)

% Extracts the differen part between two phases. 
%
% [diff_phas,diff_phas2] = isequal_phase(phases,phases2) % complete call
%
%
% INPUTS:
%
% phases: (n_phasesx5) first set of disjoints phases modelling an interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
% 
% phases2: (n2_phasesx5) second set of disjoints phases modelling an interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
%
%
% OUTPUTS:
%
% diff_phas: (n3_phasesx5) different part of the first phase. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
% 
% diff_phas2: (n4_phasesx5) different part of the second phase. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 16/Nov/06
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

% Initialization

diff_phas = phases;
diff_phas2 = phases2;


% Comparison

sp = sp1;
for i=sp(1):-1:1,
    [exists,pph2,kk,in1,in2] = evalue_phase(diff_phas2,diff_phas(i,4),diff_phas(i,5));
    if exists & isequal_phase(diff_phas(i,:),pph2),
            diff_phas = [diff_phas(1:i-1,:);diff_phas(i+1:end,:)];
            diff_phas2 = [diff_phas2(1:in1-1,:);diff_phas2(in2+1:end,:)];
    end
end
   