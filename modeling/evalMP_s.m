function [sc,e] = evalMP_s(mp_model)

% Divides the original data used to calibrate a MP model in modelled part (T) 
% and residuals (E) according to the PCA equation X = T*P'+ E.  
%
% mp_model = evalMP_s(mp_model)  % complete call
%
%
% INPUTS:
%
% mp_model: (structure) MP model (use the command "help info" for
% more info)
%
%
% OUTPUTS:
%
% sc: (KxJxI) part modelled, K(sampling times) x J(variables)
%   x I(batches)
%
% e: (KxJxI) residuals, K(sampling times) x J(variables)
%   x I(batches)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/09
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

sc = zeros(s);
e = zeros(s);

% Main code

xcs = preprocess3D(mp_model.arg.xini,mp_model.arg.prep);
for ind=1:sp(1),
    ind_i = max(1,mp_model.phases(ind,4)-mp_model.phases(ind,3));
    ind_j = mp_model.phases(ind,4) - ind_i;
    xu = unfold(xcs(ind_i:mp_model.phases(ind,5),:,:),mp_model.phases(ind,3));
    [p,t] = pcamv(xu,mp_model.phases(ind,2));
    t_sm = fold(t*p',s(3),mp_model.phases(ind,3));
    e_sm = fold(xu - t*p',s(3),mp_model.phases(ind,3));
    sc(mp_model.phases(ind,4):mp_model.phases(ind,5),:,:) = t_sm(ind_j+1:end,:,:);
    e(mp_model.phases(ind,4):mp_model.phases(ind,5),:,:) = e_sm(ind_j+1:end,:,:);       
end


    