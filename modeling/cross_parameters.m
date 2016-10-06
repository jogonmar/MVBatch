function cross=cross_parameters(leave_m,blocks_r,blocks_c,fold_m,order)

% Inicialization of parameters for 3D-crossvalidation.
%
% cross=cross_parameters % lnso cross-validation method with standart
%                          parameters
% cross=cross_parameters(leave_m,blocks_r,blocks_c,fold_m,order)  % Complete call
%
%
% INPUTS:                                                         
%
% leave_m: (text) cross-validation procedure:
%   'rkf': row-wise k-fold cross-validation.
%   'ekf': sample-wise k-fold cross-validation (by default).
%   'iekf': iterative sample-wise k-fold cross-validation.
%   'cekf': cross-corrected sample-wise k-fold cross-validation. 
%   'ckf': column-wise k-fold cross-validation (by default). 
%
% blocks_r: (1x1) maximum number of blocks of samples (Inf by default)
%
% blocks_c: (1x1) maximum number of blocks of variables (Inf by default).
%
% fold_m (text) folding method:
%   'mean': mean of all the values of a variable.
%   'first': first value (by default).
%
% order: (structure) to define a constant random ordering of columns and
%   rows.
%       order.input: (boolean)
%           true: take ordering from the structure.
%           false: compute it ramdomly (by default).
%       order.cols: (1xn_cols) columns ordering.
%       order.rows: (1xn_rows) rows ordering.
%
%
% OUTPUTS:
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
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 14/Sep/09
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

if nargin < 1, leave_m = 'ckf'; end;
if nargin < 2, blocks_r = Inf; end;
if nargin < 3, blocks_c = Inf; end;
if nargin < 4, fold_m = 'first'; end;
if nargin < 5, order.input = false; end;

% Main code

cross=struct('leave_m',leave_m,'blocks_r',blocks_r,'blocks_c',blocks_c,'fold_m',fold_m,'order',order);