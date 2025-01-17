function [obs_rec,t_est] = miss2D(x,obs,pc,miss_m,prep)

% Missing data recovery using a complete calibration data set. 
%
% [obs_rec,t_est] = miss2D(x,obs,pc,miss_m,prep)
%
% INPUTS:
%
% x: Two-way data matrix, N(observations) x M(variables)
%
% obs: uncomplete observation, 1 x M(variables)
%
% pc: number of principal components.
%
% miss_m: missing data recovery procedure:
%   'dir': direct method.
%   'iter': iterative method.
%   'tsr': trimmed score regression method.
%
% prep:
%       0: no preprocessing.
%       1: mean-centering (default).
%       2: auto-scaling.    
%
%
% OUTPUTS:
%
% obs_rec: recovered uncomplete observation, 1 x M(variables).
%
% t_est: estimated scores, 1 x pc.
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 27/Aug/09.

% Arguments checking

if nargin < 3, error('Error in the number of arguments.'); end;

if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;

if ndims(obs)~=2, error('Incorrect number of dimensions of obs.'); end;
s_o = size(obs);
if s_o(1)~=1, error('Incorrect content of obs.'); end;
if s_o(2)~=s(2), error('Incorrect content of obs.'); end;

if nargin < 4, miss_m = 'tsr'; end;
if nargin < 5, prep = 2; end;

if pc<0, error('Incorrect value of pc.'); end;


% Preprocessing

obs_rec = obs;
t_rec = [];
res_rec = [];
van = find(~isnan(obs));
vnan = find(isnan(obs));

[ccs,av,st] = preprocess2D(x,prep);
        
scs = (obs_rec-av)./st;
scs(vnan) = 0;
    
% Missing data recovery

if pc > 0, 
    
    %[p,t] = pca(ccs,pc);
    [p,t] = pca(ccs,'Centered',false,'NumComponents',pc);
    
    switch lower(miss_m)
    
        case 'dir',
            t_est = scs*p;
            scs_rec = t_est*p';

        case 'iter',
            conv = false;
            cont = 0;
            while ~conv,
                cont = cont+1;
                if cont>30, conv=true; end

                t_est = scs*p;
                scs_rec = t_est*p';
                if sum((scs(vnan)-scs_rec(vnan)).^2)<1e-10,
                    conv=true;
                else
                    scs(vnan)=scs_rec(vnan);
                end
            end
                    
        case 'tsr',
                
            TS = ccs(:,van)*p(van,:);
            B = inv(TS'*TS)*TS'*t;
            
            t_est = scs(van)*p(van,:)*B;
            scs_rec = t_est*p';
            
            if isempty(van),
                [ccs2,av2] = preprocess2D(ccs,1);
                scs_rec = va2;
            end              
            
        otherwise
               error('Incorrect miss_m.');

    end
              
else 
    scs_rec(vnan) = 0;
end

scs_rec(van) = scs(van);
obs_rec = scs_rec.*st+av;

