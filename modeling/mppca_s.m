function [mp_model,text_tot] = mppca_s(x,lag,T,absolute,gamma,minsize,n,prep,cross,console,text_tot)

% Top-Down step of MPPCA. Division in the samples of the unfolded data. 
% The improvement measure is the summed square error.
%
% [mp_model,text_tot] = mppca_s(x) % variable-wise MP model with standard parameters
% [mp_model,text_tot] = mppca_s(x,lag) % batch dynamic MP model with standard parameters
% [mp_model,text_tot] = mppca_s(x,lag,T,absolute,gamma,minsize,n,prep,cross) % output
%                                                           in MATLAB console
% [mp_model,text_tot] = mppca_s(x,lag,T,absolute,gamma,minsize,n,prep,cross,console,text_tot) % complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
% one in the row of the unfolded matrix (0 by default).
%
% T: (1x1) improvement threshold for subdivision (0 by default).
%
% absolute: (boolean) absolute improvement (true) or relative improvement
%   (false, by default).
%
% gamma: (1x1) factor to adjust the improvement of an additional PC and of a
%   division:
%       - [0-Inf]: constant value (1 by defect). 
%       - [-Inf-0): criterium of parsimony.
%
% minsize: (1x1) minimum number of sampling times in a phase (1 by default).
%
% n: (1x1) initial number of PCs (0 by default).
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: trajectory centering (average trajectory subtraction)
%       2: 1 + trajectory-scaling (scales data so that each pair variable and 
%           sampling time has variance 1) (default)  
%       3: 1 + variable-scaling (scales data so that each variable has
%           variance 1)
%       4: variable centering (subtraction of the average value of each
%           variable)
%       5: 4 + variable-scaling. 
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
% console: (1x1) handle of the EditText of the interface, 0 stands for the
%   MATLAB console (by default)
%
% text_tot: (text) input text with information of the analysis ([] by
% default).
%
%
% OUTPUTS:
%
% mp_model: (structure) MP model (use the command "help info" for
%       more info)
%
% text_tot: (text) output text with information of the analysis.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 31/Oct/08
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
if nargin < 2, lag = 0; end;
if nargin < 3, T = 0; end;
if nargin < 4, absolute = false; end;
if nargin < 5, gamma = 1; end;
if nargin < 6, minsize = 1; end;
if nargin < 7, n = 0; end;
if nargin < 8, prep = 2; end;
if nargin < 9, cross = cross_parameters; end;
if nargin < 10, console = 0; end;
if nargin < 11, text_tot = []; end;

if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if lag<0, error('Incorrect value of lag.'); end;
if lag>s(1)-1, lag=s(1)-1; end;
if (minsize<1||minsize>s(1)), error('Incorrect value of minsize.'); end;
if n<0, error('Incorrect value of n.'); end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;


% Initialization

if cross.order.input==false,
    cross.order.input=true;
    cross.order.cols=rand(1,s(2)*s(1));
    cross.order.rows=rand(1,s(3));
end

arg=struct('xini',x,'lag',lag,'T',T,'cross',cross,'absolute',absolute,'gamma',gamma,'minsize',minsize,'n',n,'prep',prep); 
pc = n;
clu = [ones(s(1),1)];
mp_model=struct('type','SW-Div','arg',arg,'pcs',n,'clu',clu);


% Greedy search

cumpress = crossval3D_s(x,n,lag,clu,cross.leave_m,cross.blocks_r,cross.blocks_c,cross.fold_m,prep,cross.order);
           
[mp_model,text_tot] = high_mppca_s(arg,clu,pc,cumpress,console,text_tot);

mp_model.phases = reduce(mp_model.tree);

mp_model = crossvalMP_s(mp_model);
