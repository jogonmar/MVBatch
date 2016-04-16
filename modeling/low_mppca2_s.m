function q2 = low_mppca2_s(x,pc,lag,minsize,prep,leave_m)

% Top-Down step, low-level procedure of MPPCA, version 2. Division in the samples
% of the unfolded data. 
%
% q2 = low_mppca2_s(x,pc) % Variable-wise data
% q2 = low_mppca2_s(x,pc,lag) % Batch dynamic data
% q2 = low_mppca2_s(x,pc,lag,minsize,prep,leave_m) % Complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%     x I(batches)
%
% pc: (1x1) number of PCs
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
%   one in the row of the unfolded matrix (0 by default).
%
% minsize: (1x1) minimum number of sampling times in a phase (1 by default).
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
% leave_m: (text) cross-validation procedure:
%   'rkf': row-wise k-fold cross-validation.
%   'skf': sample-wise k-fold cross-validation (by default).
%   'iskf': iterative sample-wise k-fold cross-validation.
%   'cskf': cross-corrected sample-wise k-fold cross-validation. 
%
%
% OUTPUTS:
%
% q2: (Kx1) Index of convenience of division in each sampling time.
%
%
% codified by: José Camacho Páez.
% version: 2.0
% last modification: 20/May/09.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, lag = 0; end;
if nargin < 4, minsize = 1; end;
if nargin < 5, prep = 2; end;
if nargin < 6, leave_m = 'skf'; end;

if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of pc.'); end;
if lag<0, error('Incorrect value of lag.'); end;
if (minsize<1||minsize>s(1)), error('Incorrect value of minsize.'); end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;

% Initialization

q2 = 0;
if s(1)<2*minsize,
    return
end

% Preprocessing

[xce,av,st] = preprocess3D(x,prep);

% Extensive search

mejor=Inf;
for i=(minsize+lag):(s(1)-minsize),   
    a =  crossval3D_s(x((i-lag-minsize+1):(i+minsize),:,:),pc,lag,[zeros(lag,1);ones(2*minsize,1)],leave_m,Inf,Inf,'first',prep);
    q2(i-lag) = (a-crossval3D_s(x((i-lag-minsize+1):(i+minsize),:,:),pc,lag,[zeros(lag,1);1*ones(minsize,1);2*ones(minsize,1)],leave_m,Inf,Inf,'first',prep))/a;
end;



