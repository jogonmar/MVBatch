function [phases2, changes] = chph_parsimony(phases,mp_group,ini,fin,Tc,absolute,val)

% Selection of the best MP sub-model for a time interval. Criterium of parsimony in the model.
%
% [phases2, changes] = chph_parsimony(phases,mp_group,ini,fin,Tc,absolute,val) % complete call
%
%
% INPUTS:
%
% phases: (n_phasesx5) original set of disjoints phases modelling the interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time].
%
% mp_group: (cell) group of MP models from which the best sub-model for the interval
%   will be obtained.
%
% ini: (1x1) initial sampling time of the interval.
%
% fin: (1x1) final sampling time of the interval.
%
% Tc: (1x1) penalization.
%
% absolute: (boolean) absolute improvement (true) or relative improvement (false).
%
% val: (1x1) value used as a baseline for absolute improvement.
%
%
% OUTPUTS:
%
% phases2: (n2_phasesx5) resulting set of disjoints phases modelling the interval. 
%       Each row contains the information of a phase, namely 
%       [PRESS, PCs, lags, initial time, end time]. 
%
% changes: (boolean) 'true' for phases2 different to phases.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Abr/07
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

% Parameters cheking

if nargin < 7, error('Error in the number of arguments.'); end;

% Initialization

changes = false;
phases2 = phases;
s = size(mp_group{1}.arg.xini);

% Search

[exists,phases_a,cumpress_a,ind1,ind2]=evalue_phase(phases2,ini,fin);
if ~exists, return; end

len=length(mp_group);        
for i=1:len,
    [exists,phases_b]=evalue_phase(mp_group{i}.phases,ini,fin);
    if exists && ~isequal_phase(phases_a,phases_b), 
        [diff_a,diff_b] = diff_phase(phases_a,phases_b);       
        npa=number_params(diff_a,s(2));
        npb=number_params(diff_b,s(2));
        cumpress_da = sum(diff_a(:,1));
        cumpress_db = sum(diff_b(:,1));
        
        T = Tc*(npb-npa)/max(npb,npa);
        
        if absolute,
            sda=size(diff_a);
            baseline=0;
            for i=1:sda(1),
                baseline = baseline + sum(val(diff_a(i,4):diff_a(i,5)));
            end
        else
            baseline=max(cumpress_da,cumpress_db);
        end
                
        if (cumpress_da-cumpress_db)/baseline> T,
            phases_a=phases_b;
            changes=true;
        end
    end
end    

phases2=[phases2(1:ind1-1,:);phases_a;phases2(ind2+1:end,:)];
