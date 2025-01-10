function [p,t] = pcamv(x,pc)

% Principal Component Analysis.
%
% [p,t] = pca(x,pc)     % complete call
%
% INPUTS:
%
% x: (NxM) Two-way batch data matrix, N(observations) x M(variables)
%
% pc: number of principal components.
%
%
% OUTPUTS:
%
% p: (M x pc) matrix of loadings.
%
% t: (N x pc) matrix of scores.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 10/Jan/2025
%
% Copyright (C) 2025  University of Granada, Granada
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

%[p,t] = pca_pp(x,pc);
[p,t] = pca(x,'Centered',false,'NumComponents',pc);

p = p(:,end);
t = t(:,end);


