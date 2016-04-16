function [clu,q1,q2,indc] = low_mppca_s(x,pc,lag,minsize,prep)

% Top-Down step, low-level procedure of MPPCA. Division in the samples
% of the unfolded data. The improvement measure is the summed square error.
%
% [clu,q1,q2,indc] = low_mppca_s(x,pc) % Variable-wise data
% [clu,q1,q2,indc] = low_mppca_s(x,pc,lag) % Batch dynamic data
% [clu,q1,q2,indc] = low_mppca_s(x,pc,lag,minsize,prep) % Complete call
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
% minsize: (1x1) minimum number of sampling times in a phase (1 by
%   default).
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
%
% OUTPUTS:
%
% clu: (K-lagx1) vector with the assignment of the sampling times to the phases.
%
% q1: (1x1) SSE of complete model.
%
% q2: (Kx1) SSE of divided model.
%
% indc: (1x1) last sampling time of first phase.
%
%
% codified by: José Camacho Páez.
% version: 1.2
% last modification: 20/Aug/09.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, lag = 0; end;
if nargin < 4, minsize = 1; end;
if nargin < 5, prep = 2; end;

if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of pc.'); end;
if lag<0, error('Incorrect value of lag.'); end;
if (minsize<1||minsize>s(1)), error('Incorrect value of minsize.'); end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;

% Initialization

clu = ones(s(1)-lag,1);
q1 = 0;
q2 = 0;
indc = 0; 
if s(1)<2*minsize,
    return
end

% Preprocessing

[xce,av,st] = preprocess3D(x,prep);

% Complet model

xu=unfold(xce,lag);
[p,t] = pca(xu,pc);
q1 = sum(sum((xu-t*p').^2));

% Extensive search
mejor=Inf;
for i=minsize:(s(1)-lag-minsize),    
    
    [pa,ta] = pca(xu(1:i*s(3),:),pc);
    [pb,tb] = pca(xu(i*s(3)+1:(s(1)-lag)*s(3),:),pc);
    
    q2(i+lag) = sum(sum((xu(1:i*s(3),:)-ta*pa').^2)) + sum(sum((xu(i*s(3)+1:(s(1)-lag)*s(3),:)-tb*pb').^2));
    
    if q2(i+lag) < mejor,
        mejor = q2(i+lag);
        indc = i;
    end;
end;

clu = [ones(indc,1);2*ones(s(1)-lag-indc,1)];

