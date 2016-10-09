function [EvolvingDcontribution,EvolvingQcontribution]=EvolvingContributionToDQ(xini,test,phases,prep)

% Compute the evolving contribution to the D and Q statistics for a test batch from the beginning to the end of the run. 
% The procedure is based on the original proposal publised in Technometrics, 1995, 37:1.
%
% CALL:
%
% [EvolvingDcontribution,EvolvingQcontribution]=EvolvingContributionToDQ(xini,test,phases,prep) % complete call
%
% INPUTS:
%
% xini: (KxJxI) three-way batch data array for calibration, K(sampling times) 
%       x J(variables) x I(batches)
%
% test: (KxJ) two-way batch data matrix for test, K(sampling times) 
%       x J(variables) 
%
% phases: (n_phasesx5) phases of the MP model. Each row contains the information 
%   of a phase, namely [PRESS, PCs, lags, initial time, end time]. 
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
% OUTPUTS:
%
% EvolvingDcontribution: (KxJ) contribution to the D statistic for all the
%                        time points when monitoring the batch in an evolving
%                        fashion. i.e. the variance-covariance matrix is 
%                        estimated from the start of the batch to the current 
%                        time point k. 
%
% EvolvingQcontribution: (KxJ) contribution to the Q statistic for all the
%                        time points when monitoring the batch in an evolving
%                        fashion. i.e. the variance-covariance matrix is 
%                        estimated from the start of the batch to the current 
%                        time point k.
%
% Copyright (C) 2016  Technical University of Valencia, Valencia
% Copyright (C) 2016  José M. González Martínez
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


%% Arguments checking
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

if ndims(xini)~=3, error('Incorrect number of dimensions of xini.'); end;
s = size(xini);
if find(s<1), error('Incorrect content of xini.'); end;
if ndims(test)~=2, error('Incorrect number of dimensions of test.'); end;
st=size(test);
if s(1)~=st(1) || s(2)~=st(2),
    error('Incorrect content of test.')
end
if ndims(phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp=size(phases);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of phases.'); end;
if find(phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(phases(:,4:5)<1), error('Incorrect content of phases.'); end;
if find(phases(:,3:5)>s(1)), error('Incorrect content of phases.'); end;
if nargin < 4, prep = 2; end;

%% Initialize parameters
EvolvingDcontribution = [];
EvolvingQcontribution = [];


%% Main code
m=sp(1);

% Preprocess the calibration data set
[xce,av,sta] = preprocess3D(xini,prep);

% Preprocess the test batch with the calibration preprocessing parameters
teste = (test-av)./sta;

if phases(:,2)>0,
    for i=1:m,
        ind=(max(phases(i,4)-phases(i,3),1):phases(i,5));
        ind_ini=find(ind==phases(i,4));
        xu=unfold(xce(ind,:,:),phases(i,3));
        [U,S,V] = svd(xu,'econ');
        tAll = U*S;
        pAll = V;
        p = pAll(:,1:phases(i,2));
        t = tAll(:,1:phases(i,2));
        
        % Projection of the test batch onto the latent subspace
        testu=unfold(teste(ind,:,:),phases(i,3));
        tpred = testu*p;
        resb=testu-tpred*p';
        resb=fold(resb,1,phases(i,3));       
        if ind_ini<phases(i,3)+1,
            menor_en=phases(i,3)-ind_ini;
            % Estimate covariance matrices for batch trajectory forecasting
            % with TDR
            theta = cov(tAll);
            theta_A = cov(t);

            for j=0:menor_en,
                jindb=1:s(2)*(ind_ini+j);
                jind2=phases(i,4)+j; 
                % Trajectory imputation with TSR
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*xu(1:s(3),jindb)';
                cov_inv=inv(cov(t_t'));
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*testu(1,jindb)';
                % Calculate contributions to the D statistic
                ck=[];
                for z = jindb(end-s(2)+1):jindb(end)
                    dist=0;
                    for a=1:phases(i,2)
                        dist =  dist  + ((t_t(a)/sqrt(cov_inv(a,a)))^2)*p(z,a)^2;
                    end 
                    ck = [ck, testu(1,z).*sqrt(dist)];

                end
                EvolvingDcontribution = [EvolvingDcontribution; ck];
                EvolvingQcontribution = [EvolvingQcontribution; (permute(teste(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-s(2)+1:end),:)').^2]; 
            end        
        end
        ssc=size(t);
        j=1;
        % Section to estimate the statistics for variable-wise models and
        % for batch-wise/dynamic-wise in their K-lags samples left
        for o=1:s(3):ssc(1),
            sc_model = t(o:o+s(3)-1,:);
            cov_inv = inv(cov(sc_model));
            indb=size(testu,2)-s(2)+1;
            ck=[];
            for z=indb:size(testu,2)
                dist=0;
                for a=1:phases(i,2)
                    dist =  dist  + ((tpred(j,a)/sqrt(cov_inv(a,a)))^2)*p(z,a)^2;
                end 
                ck = [ck, testu(j,z).*sqrt(dist)];
            end
            EvolvingDcontribution = [EvolvingDcontribution; ck];            
            j=j+1;
       end
        EvolvingQcontribution = [EvolvingQcontribution; resb((phases(i,3)+1):end,:,:)]; 
    end
end






