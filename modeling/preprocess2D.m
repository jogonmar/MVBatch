function [xp,average,scale] = preprocess2D(x,prep)

% Preprocess 2-way data.
%
% [xp,average,scale] = preprocess2D(x)          % for mean centering
% [xp,average,scale] = preprocess2D(x,prep)     % complete call
%
% INPUTS:
%
% x: (NxM) Two-way batch data matrix, N(observations) x M(variables)
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)  
%
% OUTPUTS:
%
% xp: (NxM) preprocessed data.
%
% average: (1 x M) sample average according to the preprocessing method.
%
% scale: (1 x M) sample scale according to the preprocessing method.
%
%
% codified by: José Camacho Páez.
% last modification: 29/Jul/09.

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 2, prep = 1; end;
if (prep<0||prep>2), error('Incorrect value of prep.'); end;

% Computation

switch prep,
    
    case 1, % mean centering
        
        for i=1:size(x,2),
            v = find(~isnan(x(:,i)));
            average(i) = sum(x(v,i),1)/length(v);
        end
        scale = ones(1,s(2));
        xp = x - ones(s(1),1)*average;
        
    case 2, % Trajectory centering and scaling
        
        for i=1:size(x,2),
            v = find(~isnan(x(:,i)));
            average(i) = sum(x(v,i),1)/length(v);
            xc(:,i) = x(:,i) - ones(s(1),1)*average(i); 
            scale(i) = sqrt(sum(xc(v,i).^2,1)/(length(v)-1));
        end
        scale(find(scale==0))=1; 
        xp = xc./(ones(s(1),1)*scale);
        
    otherwise, % No preprocessing 
        average = zeros(1,s(2));     
        scale = ones(1,s(2)); 
        xp = x;
end


