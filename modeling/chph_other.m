function [phases2, changes] = chph_other(phases,mp_group,ini,fin,T_s,T_lv,T_lmv,absolute,val)

% Selection of the best sub-model for an interval. Several criteria.
%
% [phases2, changes] =
% chph_other(phases,mp_group,ini,fin,T_s,T_lv,T_lmv,absolute,val)
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
% T_s: (1x1) improvement threshold for subdivisions.
%
% T_lv: (1x1) improvement threshold for including latent variables.
%
% T_lmv: (1x1) improvement threshold for including lagged measurement-vectors.
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
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 19/Abr/07.

% Parameters checking

if nargin < 9, error('Error in the number of arguments.'); end;

% Initialization

changes = false;
phases2 = phases;

% Search

[exists,phases_a,cumpress_a,ind1,ind2]=evalue_phase(phases2,ini,fin);
if ~exists, return; end

len=length(mp_group);        
for i=1:len,
    [exists,phases_b]=evalue_phase(mp_group{i}.phases,ini,fin);
    if exists && ~isequal_phase(phases_a,phases_b), 
        [diff_a,diff_b] = diff_phase(phases_a,phases_b);       
        spa=size(diff_a);
        spb=size(diff_b);
        cumpress_da = sum(diff_a(:,1));
        cumpress_db = sum(diff_b(:,1));
        np=spb(1)-spa(1); % difference in the number of phases
        npc=sum(diff_b(:,2))/spb(1)-sum(diff_a(:,2))/spa(1); % difference in the averaged PCs
        nl=sum(diff_b(:,3))/spb(1)-sum(diff_a(:,3))/spa(1); % difference in the averaged LMV
        
        if isnan(nl),
            nl=0;
        end
        
        T = T_s*np+T_lv*npc+T_lmv*nl;
        
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
