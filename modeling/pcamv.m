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
% codified by: José Camacho Páez.
% last modification: 23/Apr/09.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of prep.'); end;
dmin = min(s);
if pc>dmin, pc=dmin; end;

% Computation

if 10*s(1)>s(2),
        [p,t]=princomp(x);
        p = p(:,1:pc);
        t = t(:,1:pc);
else,
        [p,t]=princomp(x,'econ');
        p = p(:,1:pc);
        t = t(:,1:pc);
end
        



