function mp_model = crossvalMP_s(mp_model,cross)

% Cross-validation for square-prediction-errors computing in Multi-Phase
% models of 3-to-2 way batch data, divisions in the samples.  
%
% mp_model = crossvalMP_s(mp_model)  % cross-validation parameters in the
%                                       model structure.
% mp_model = crossvalMP_s(mp_model,cross)  % complete call
%
%
% INPUTS:
%
% mp_model: (structure) MP model (use the command "help info" for more info)
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
% mp_model: (structure) MP model with the values of CUMPRESS and PRESS.
%   (use the command "help info" for more info)
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 12/Abr/07
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
if nargin < 2, cross = mp_model.arg.cross; end;

if ndims(mp_model.arg.xini)~=3, error('Incorrect number of dimensions of xini.'); end;
s = size(mp_model.arg.xini);
if find(s<1), error('Incorrect content of xini.'); end;
if (mp_model.arg.prep<0||mp_model.arg.prep>5), error('Incorrect value of prep.'); end;

if ndims(mp_model.phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp = size(mp_model.phases);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,4:5)<1), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,3:5)>s(1)), error('Incorrect content of phases.'); end;


% Initialization

mp_model.pem = zeros(s);
mp_model.press = zeros(s(1),1);
mp_model.cumpress = 0;

% Cross-validation

for ind=1:sp(1),
    
    clu=zeros(s(1),1);
    clu(mp_model.phases(ind,4):mp_model.phases(ind,5)) = 1;
    [cumpress_inc,press,pem] = crossval3D_s(mp_model.arg.xini,mp_model.phases(ind,2),mp_model.phases(ind,3),clu,cross.leave_m,cross.blocks_r,cross.blocks_c,cross.fold_m,mp_model.arg.prep,cross.order);    
    mp_model.pem(mp_model.phases(ind,4):mp_model.phases(ind,5),:,:) = pem(mp_model.phases(ind,4):mp_model.phases(ind,5),:,:);
             
    if isequal(lower(mp_model.arg.cross.fold_m),'mean') && (mp_model.phases(ind,4)>mp_model.phases(ind,3)+1),
        ind_p=(mp_model.phases(ind,4)-mp_model.phases(ind,3):mp_model.phases(ind,4)-1);
        w1t=(length(ind_p):-1:1)'*ones(1,s(2));        
        w1=min(w1t,mp_model.phases(ind-1,3)+1);
        
        w2=(1:length(ind_p))'*ones(1,s(2));
        w2=min(w2,length(mp_model.phases(ind,4):mp_model.phases(ind,5)));
            
        for j=1:s(3); 
            mp_model.pem(ind_p,:,j) = (w1.*mp_model.pem(ind_p,:,j) + (w1t+w2).*pem(ind_p,:,j))./(w1+w2);      
        end      
    end        
    
end

mp_model.press = sum(sum(mp_model.pem.^2,3),2);

for ind=1:sp(1),
    mp_model.phases(ind,1) = sum(mp_model.press(mp_model.phases(ind,4):mp_model.phases(ind,5)));
end

mp_model.cumpress = sum(mp_model.press);
    