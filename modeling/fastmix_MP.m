function [best,mp_group2] = fastmix_MP(mp_group,Tm,criterium,cross_group,cross)

% Heuristic algorithm for the mixture of several MPPCA models. 
%
% [best,mp_group2] = fastmix_MP(mp_group) % standard parameters
% [best,mp_group2] = fastmix_MP(mp_group,Tm,criterium,cross_group,cross) % complete
%                                                                          call
%
% INPUTS:
%
% mp_group: (cell) containing a set of MPPCA models.
%
% Tm: (1x1) improvement threshold for the addition of components or sub-models (0 by default).
%   When criterium is 'all', Tm should be a (3x1) vector.
%
% criterium: (text) criterium for phases comparison:
%       - 'parsimony': achieve a parsimonious MP model.
%           Threshold = Tm * (n_parB - n_parA)/max(n_parB - n_parA).
%       - 'parsimony2': achieve a parsimonious covariance matrix.
%           Threshold = Tm * (n_par_covB - n_par_covA)/max(n_par_covB - n_par_covA).
%       - 'LMV': reduce the number of Lagged Measurement Vectors.
%           Threshold = Tm * (average_LMVB - average_LMVA).
%       - 'phases': reduce the number of phases.
%           Threshold = Tm * (n_phaB - n_phaA).  
%       - 'pcs': reduce the number of principal components.
%           Threshold = Tm * (average_pcsB - average_pcsA). 
%       - 'all': reduce the number of LMVs, phases and principal components.
%           Threshold = Tm(1) * (average_LMVB - average_LMVA) + Tm(2) * (n_phaB - n_phaA)
%                       + Tm(3) * (average_pcsB - average_pcsA). 
%
% cross_group: (boolean) use cross-validation for comparisons (true by default).
%
% cross: (structure) the parameters.  
%   cross.leave_m: (text) cross-validation procedure
%   cross.blocks_r: (1x1) maximum number of blocks of samples
%   cross.blocks_c: (1x1) maximum number of blocks of variables
%   cross.fold_m: (text) folding method
%   cross.order: (structure) to define a constant random ordering of columns and
%       rows.
%       cross.order.input: (boolean)
%           true: take ordering from the structure.
%           false: compute it ramdomly (by default).
%       cross.order.cols: (1xn_cols) columns ordering.
%       cross.order.rows: (1xn_rows) rows ordering.
%
%
% OUTPUTS:
%
% best: (structure) MP model (use the command "help info" for
%       more info) resulting from the best mixture of MPPCA models.
%
% mp_group2: (cell) total set of MPPCA models generated.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/Feb/09
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
if nargin < 2, Tm=0; end;
if nargin < 3, criterium='parsimony'; end;
if nargin < 4, cross_group=true; end;

if length(mp_group)<1, error('Incorrect content of mp_group.'); end;

if isequal(criterium,'all'),
    if length(Tm)~=3, 
        error('Incorrect length of Tm.');
    end
else
    if length(Tm)~=1, 
        error('Incorrect length of Tm.');
    end
end

if nargin < 5, cross = mp_group{1}.arg.cross; end;


% Initialization

arg = mp_group{1}.arg;
arg.Tm = Tm;
type = mp_group{1}.type;

mp_group2 = mp_group;

for i=1:length(mp_group), % Obtaining the phases with the new threshold of division 
    mp3=thres_subTD(mp_group{i},Tm);    
    already=false;
    for j=1:length(mp_group2), % Only add if it was not yet
        if isequal_phase(mp_group2{j}.phases,mp3.phases),
            already=true;
            j=length(mp_group2)+1;
        end
    end
    if ~already,
        mp_group2{end+1}=mp3;
    end
    
    mp4=thres_subDT(mp_group{i},Tm);
    already=false;
    for j=1:length(mp_group2), % Only add if it was not yet
        if isequal_phase(mp_group2{j}.phases,mp4.phases),
            already=true;
            j=length(mp_group2)+1;
        end
    end
    if ~already,
        mp_group2{end+1}=mp4;
    end
end

s=size(arg.xini);   
[cumpress_0,press_0] = crossval3D_s(arg.xini,0,0,ones(s(1),1),cross.leave_m,cross.blocks_r,cross.blocks_c,cross.fold_m,arg.prep,cross.order);

% Mixture

changes3=true;
len_prev=0;
while changes3, % while adding new mixtures
  
    len=length(mp_group2);
    total=len;
    changes3=false;
    
    for u=len_prev+1:len, % current models used to initialize mixtures
         
        sp=size(mp_group2{u}.phases); % take model M
        phases=mp_group2{u}.phases;
 
        changes=false;
        for i=1:sp(1),
            for o=1:sp(1)-i+1,
                ini=mp_group2{u}.phases(o,4);
                fin=mp_group2{u}.phases(o+i-1,5);  
                
                switch lower(criterium)
                    
                    case 'parsimony',
                        [phases,changes2]=chph_parsimony(phases,mp_group2,ini,fin,Tm,arg.absolute,press_0);
                    case 'parsimony2',
                        [phases,changes2]=chph_parsimony2(phases,mp_group2,ini,fin,Tm,arg.absolute,press_0);
                    case 'phases',
                        [phases,changes2]=chph_other(phases,mp_group2,ini,fin,Tm,0,0,arg.absolute,press_0); 
                    case 'lmv',
                        [phases,changes2]=chph_other(phases,mp_group2,ini,fin,0,0,Tm,arg.absolute,press_0); 
                    case 'pcs',
                        [phases,changes2]=chph_other(phases,mp_group2,ini,fin,0,Tm,0,arg.absolute,press_0);  
                    case 'all',
                        [phases,changes2]=chph_other(phases,mp_group2,ini,fin,Tm(1),Tm(2),Tm(3),arg.absolute,press_0);   
                    otherwise
                        error('Incorrect criterium.');
                end
                if changes2,
                    changes=true;
                end
            end
        end

        total=length(mp_group2);
        if changes,
                                                       
           already = false; 
           for j=1:total,
               if isequal_phase(mp_group2{j}.phases,phases),
                   already = true;
               end
           end
        
           if ~already, 
               [phases,clu]=reduce(phases);
               mp_model=struct('type',type,'arg',arg,'clu',clu,'phases',phases,'tree',phases);  
                        
               if cross_group,
                   if isequal(type,'SW-Div'),
                       mp_model = crossvalMP_s(mp_model,cross);
                   else
                       mp_model = crossvalMP_v(mp_model,cross);
                   end
               end

               mp_group2{total+1}=mp_model;
               changes3=true;
                    
            end
        end 
        best = phases;
    end
    len_prev = len;
    len = total;
end

for i=length(mp_group2):-1:1,
    if isequal_phase(mp_group2{i}.phases,best),
        if ~cross_group,
            if isequal(type,'SW-Div'),
                best = crossvalMP_s(mp_group2{i},cross);
            else
                best = crossvalMP_v(mp_group2{i},cross);
            end
        else
            best = mp_group2{i};
        end        
        best.arg.Tm=Tm;
        return;
    end
end

            
