function [evolD,evolQ,D,Q,EvolvingDcontribution,EvolvingQcontribution] = multiphaseProjection(xini,test,phases,P,T,av,stnd,flagCont)

% Computes the D-statistic and Q-statistics values for on-line monitoring.
%
% [evolD,evolQ,D,Q,res,EvolvingDcontribution,EvolvingQcontribution] = multiphaseProjection(xini,test,phases,P,T,av,stnd) % call
% [evolD,evolQ,D,Q,res,EvolvingDcontribution,EvolvingQcontribution] = multiphaseProjection(xini,test,phases,P,T,av,stnd) % complete call
%
% INPUTS:
%
% xini: (KxJxI) three-way batch data matrix for calibration, K(sampling times) 
%       x J(variables) x I(batches)
%
% test: (KxJ) three-way batch data matrix for test, K(sampling times) 
%       x J(variables)
%
% phases: (n_phasesx5) phases of the MP model. Each row contains the information 
%   of a phase, namely [PRESS, PCs, lags, initial time, end time]. 
%
% P: {phases x 1} cell array of [J_phase x A_phase] matrices of loadings for all phases.
%
% T: {phases x 1}cell array of [I_phase x A_phase] matrices of scores for all phases. 
%
% mn: [KxJ] matrix of averages.
%
% stnd: [KxJ] matrix of standard deviations.
%
% flagCont: (1x1) boolean to indicate whether the contributions to D and Q
% statistics must be computed.
%
% OUTPUTS
%
% evolD: (Kx1) online D-statistic values for the test batch.
%
% evolQ: (Kx1) online Q-statistic values for the test batch.
%
% D: (Ix1) overall D-statistic values for the test batch (if batch-wise modeling).
%
% Q: (Ix1) overall Q-statistic values for the test batch (if batch-wise modeling).
%
% coded by: José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2017  José M. González Martínez
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
assert(nargin >= 7, 'Incorrect number of input parameters');
st=size(test);
s=size(xini);
assert(length(st)<3,'Projection is done batch to batch. Please, provide a two-way array for a single batch.');
st(3)=1;
assert(s(1)==st(1) && s(2)==st(2),'Incorrect content of test');
sscores = size(T);
sloads  = size(P);
smn     = size(av);
sstd    = size(stnd);
assert(sscores(2)== sloads(2),'Incorrect number of PCs');
assert(smn(2)== sstd(2),'Incorrect number of pre-processing parameters');
assert(smn(2)== s(2),'Incorrect number of pre-processing parameters');
assert(smn(1)== s(1),'Incorrect number of averages');
assert(sstd(1)== s(1),'Incorrect number of standard deviations');
if nargin < 8, flagCont = 0; end
assert(flagCont == 0 || flagCont == 1,'Incorrect flag value to compute contribution. It must be either 0 (no compute) or 1 (compute).'); 

sp = size(phases);

% Initialization
xce = nan(s(1),s(2),s(3));
teste = nan(st(1),st(2),st(3));
evolD = nan(st(1),1);
evolQ = nan(st(1),1);
EvolvingDcontribution = nan(st(1),st(2));
EvolvingQcontribution = nan(st(1),st(2));
D = nan(st(3),1);
Q = nan(st(3),1);
res=zeros(st(3),st(2),st(1));
itersamples = 1;

% Main code
if ndims(test)==3
    for i=1:st(3)
        teste(:,:,i) = (test(:,:,i)-av)./stnd;
    end
else
    teste = (test-av)./stnd;
end

for i=1:s(3)
    xce(:,:,i) = (xini(:,:,i)-av)./stnd;
end

 for i=1:sp(1)   
    ind=(max(phases(i,4)-phases(i,3),1):phases(i,5));
    ind_ini=find(ind==phases(i,4));
    xu=unfold(xce(ind,:,:),phases(i,3));

    p = P{i}(:,1:phases(i,2));
    t = T{i}(:,1:phases(i,2));
    testu=unfold(teste(ind,:,:),phases(i,3));
    tpred = testu*p;

    resb=testu-tpred*p';
    resb=fold(resb,1,phases(i,3));
    res(:,:,ind(ind_ini:end))=permute(resb(ind_ini:end,:,:),[3 2 1]);

     %% Code for post-batch off-line process monitoring
    if phases(i,3) == st(1)-1
        Sinv = inv(cov(t));
        D = zeros(st(3),1);
        for j=1:st(3)
            D(j) = (tpred(j,:))*Sinv*(tpred(j,:))';
        end
        Q =  sum(sum(resb.^2));
    end

    if ind_ini<phases(i,3)+1
        menor_en=phases(i,3)-ind_ini;
        % Estimate covariance matrices for TSR-based imputation
        theta = cov(T{i});
        theta_A = cov(t);
        for j=0:menor_en
            jindb=1:st(2)*(ind_ini+j);
            jind2=phases(i,4)+j; 
            % TSR-based imputation
            t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*P{i,1}(jindb,:)*theta*P{i,1}(jindb,:)'*p(jindb,:))*p(jindb,:)'*xu(1:s(3),jindb)';
            cov_inv=inv(cov(t_t'));
            t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*P{i,1}(jindb,:)*theta*P{i,1}(jindb,:)'*p(jindb,:))*p(jindb,:)'*testu(1,jindb)';
            % Calculation of the D and Q statistics
            evolD(itersamples) = t_t'*cov_inv*t_t;
            evolQ(itersamples) = sum((permute(teste(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-st(2)+1:end),:)').^2);
            
            if flagCont
                % Calculate contributions to D statistic
                contrD = zeros(1,s(2));
                for z = jindb(end-s(2)+1):jindb(end)
                    dist=0;
                    for a=1:phases(i,2)
                        dist =  dist  + ((t_t(a)/sqrt(cov_inv(a,a)))^2)*p(z,a)^2;
                    end 
                    contrD(z-jindb(end-s(2)+1)+1) = testu(1,z).*sqrt(dist);
                end
                % Add the calculated contributions to D statistic to the
                % corresponding array
                EvolvingDcontribution(itersamples,:) = contrD;
                % Calculate contributions to Q statistic
                EvolvingQcontribution(itersamples,:) = (permute(teste(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-s(2)+1:end),:)').^2; 
            end
            itersamples = itersamples + 1;
        end        
    end
    ssc=size(t);
    z=1;
    
    % Compute the Q statistic for the remaining samples in batch-dynamic models
    evolQ(itersamples:phases(i,5)) = squeeze(sum(sum(resb((phases(i,3)+1):end,:,:).^2,3),2));
    % Compute contributions to Q statistic for the remaining in batch-dynamic models
    if flagCont,EvolvingQcontribution(itersamples:phases(i,5),:) = resb((phases(i,3)+1):end,:,:); end
    
    % Compute the D statistic and its contributions for the remaining samples in batch-dynamic models
    indb=size(testu,2)-s(2)+1;
    for u=1:s(3):ssc(1)
        sc_model = t(u:u+s(3)-1,:);
        cov_inv = inv(cov(sc_model));
        evolD(itersamples) = tpred(z,:)*cov_inv*tpred(z,:)';
        if flagCont
            % Compute contribution the D-statistic
            contrD = zeros(1,numel(indb:size(testu,2)),1);
            for k=indb:size(testu,2)
                dist=0;
                for a=1:phases(i,2)
                    dist =  dist  + ((tpred(z,a)/sqrt(cov_inv(a,a)))^2)*p(k,a)^2;
                end 
                contrD(k-indb+1) = testu(z,k).*sqrt(dist);
            end
            if flagCont,EvolvingDcontribution(itersamples,:) = contrD;end
        end
        z=z+1;
        itersamples=itersamples+1;
    end
    
 end      
