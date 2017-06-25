function [pdet,ph,pr,detect,h,r,ph95,pr95,h95,r95,t,q]=plot_onstat(xini, test, phases, prep, opt, alph, alpr, alph95, alpr95, s_sBn, axes1, axes2, alpoh, alpor, alpoh95, alpor95, axes3, axes4, posTest)

% Computes the D-statistic and SPE charts for on-line monitoring. 
%
% [pdet,ph,pr,detect,h,r,ph95,pr95,h95,r95,t,q] = plot_onstat(xini, test, phases, 
%   prep, opt) % call with standard parameters
%
% [pdet,ph,pr,detect,h,r,ph95,pr95,h95,r95,t,q] = plot_onstat(xini, test, phases, 
%   prep, opt, alph, alpr, alph95, alpr95) % output in MATLAB console
%
% [pdet,ph,pr,detect,h,r,ph95,pr95,h95,r95,t,q] = plot_onstat(xini, test, phases, 
%   prep, opt, alph, alpr, alph95, alpr95, axes1, axes2, alpoh, alpor, alpoh95, 
%   alpor95, axes3, axes4, posTest) % complete call
%
%
% INPUTS:
%
% xini: (KxJxI) three-way batch data matrix for calibration, K(sampling times) 
%       x J(variables) x I(batches)
%
% test: (KxJxI2) three-way batch data matrix for test, K(sampling times) 
%       x J(variables) x I2(batches)
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
% opt: boolean (1x1) 
%       true: plot results.
%       false: do not plot results.
%
% alph: imposed significance level (alpha) for the 99% confidence limit in the 
%   D-statistic. 
%
% alpr: imposed significance level (alpha) for the 99% confidence limit in the 
%   SPE.
%
% alph95: imposed significance level (alpha) for the 95% confidence limit in the 
%   D-statistic. 
%
% alpr95: imposed significance level (alpha) for the 95% confidence limit in the 
%   SPE.
%
% s_sBn: number of synchronized samples to be monitored (default: Inf, i.e. all)
%
% axes1: handle to the axes where the D-statistic chart is plotted.
%
% axes2: handle to the axes where the SPE chart is plotted.
%
% alpoh: imposed significance level (alpha) for the 99% confidence limit in the 
%   overall D-statistic. 
%
% alpor: imposed significance level (alpha) for the 99% confidence limit in the 
%   overall SPE.
%
% alpoh95: imposed significance level (alpha) for the 95% confidence limit in the 
%   overall D-statistic. 
%
% alppr95: imposed significance level (alpha) for the 95% confidence limit in the 
%   overall SPE.
%
% axes3: handle to the axes where the overall D-statistic chart is plotted.
%
% axes4: handle to the axes where the overall SPE chart is plotted.
%
% postTest: position on the graph of the values estimated for the test batch. 
%
%
% OUTPUTS:
%
% pdet: (1x1) percentage of sampling times where either the D-statistic
%   or the SPE for the test batch exceeds the corresponding 99% control 
%   limit.
%
% ph: (1x1) percentage of sampling times where the D-statistic for the test
%   batch exceeds the 99% control limit.
%
% pr: (1x1) percentage of sampling times where the SPE for the test batch
%   exceeds the 99% control limit.
%
% detect: (1xK) vector with 1s in the sampling times where either the 
%   D-statistic or the SPE for the test batch exceeds the corresponding 99%
%   control limit and 0s in the rest.
%
% h: (1xK) vector with 1s in the sampling times where the D-statistic for 
%   the test batch exceeds the 99% control limit and 0s in the rest.
%
% r: (1xK) vector with 1s in the sampling times where the SPE for the test
%   batch exceeds the 99% control limit and 0s in the rest.
%
% ph95: (1x1) percentage of sampling times where the D-statistic for the test
%   batch exceeds the 95% control limit.
%
% pr95: (1x1) percentage of sampling times where the SPE for the test batch
%   exceeds the 95% control limit.
%
% h95: (1xK) vector with 1s in the sampling times where the D-statistic for 
%   the test batch exceeds the 95% control limit and 0s in the rest.
%
% r95: (1xK) vector with 1s in the sampling times where the SPE for the test
%   batch exceeds the 95% control limit and 0s in the rest.
%
% t: (1xK) D-statistic values for the test batch.
%
% q: (1xK) SPE values for the test batch.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: 13/Dic/11
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez, José M. González Martínez
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

if nargin < 3, error('Numero de argumentos erroneos.'); end;


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
if nargin < 5, opt = 1; end;

if (prep<0||prep>5), error('Incorrect value of prep.'); end;

if nargin < 6, alph = 0.01; end;
if nargin < 7, alpr = 0.01; end;
if nargin < 8, alph95 = 0.05; end;
if nargin < 9, alpr95 = 0.05; end;
if nargin < 10, s_sBn = Inf; end
if nargin < 11 || isempty(axes1) && opt, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 12 || isempty(axes2) && opt, 
    h2 = figure;
    axes2 = axes; 
end;

if nargin < 13, alpoh = 0.01; end;
if nargin < 14, alpor = 0.01; end;
if nargin < 15, alpoh95 = 0.05; end;
if nargin < 16, alpor95 = 0.05; end;

if (nargin < 17 || isempty(axes3)) && phases(1,3) == s(1)-1 && opt, 
    h3 = figure;
    axes3 = axes; 
end;
if (nargin < 18 || isempty(axes4)) && phases(1,3) == s(1)-1 && opt, 
    h4 = figure;
    axes4 = axes; 
end;

if nargin < 19, posTest = s(3);end

% Main code

num = s(3); 
m=sp(1);

[xce,av,sta] = preprocess3D(xini,prep);

if ndims(test)==3,
    for i=1:st(3),
        teste(:,:,i) = (test(:,:,i)-av)./sta;
    end
else
    teste = (test-av)./sta;
    st(3)=1;
end

t2=[];
q=[];
pcs=[];
res=zeros(s(3),s(2),s(1));
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
        testu=unfold(teste(ind,:,:),phases(i,3));
        tpred = testu*p;
        
        resb=xu-t*p';
        ssqres = sum(resb.^2,2);
        resb=fold(resb,s(3),phases(i,3));
        res(:,:,ind(ind_ini:end))=permute(resb(ind_ini:end,:,:),[3 2 1]);


        
         %% Code for post-batch off-line process monitoring
        if phases(i,3) == s(1)-1
            Sinv = inv(cov(t));
            T2v = zeros(s(3),1);
            for j=1:s(3)
                T2v(j) = (t(j,:))*Sinv*(t(j,:))';
            end
            T2vpred = tpred *Sinv* tpred';
        end

        resb=testu-tpred*p';
        ssqrespred =  sum(resb.^2,2);
        resb=fold(resb,1,phases(i,3));
        
        if ind_ini<phases(i,3)+1,
            menor_en=phases(i,3)-ind_ini;
            % Estimate covariance matrices for TSR-based imputation
            theta = cov(tAll);
            theta_A = cov(t);
            for j=0:menor_en,
                jindb=1:s(2)*(ind_ini+j);
                jind2=phases(i,4)+j; 
                % IMPUTATION USING TSR
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*xu(1:(s(3)),jindb)';
                cov_inv=inv(cov(t_t'));
                res(:,:,jind2)=permute(xce(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-s(2)+1:end),:)';
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*testu(1,jindb)';
                % CALCULATE STATISTICS
                t2=[t2;t_t'*cov_inv*t_t];
                q=[q ;sum((permute(teste(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-s(2)+1:end),:)').^2)];
                   
            end        
        end
        ssc=size(t);
        j=1;
        for u=1:s(3):ssc(1),
            sc_model = t(u:u+s(3)-1,:);
            cov_inv = inv(cov(sc_model));
            t2=[t2;tpred(j,:)*cov_inv*tpred(j,:)'];
            j=j+1;
        end
    
        q=[q;squeeze(sum(sum(resb((phases(i,3)+1):end,:,:).^2,3),2))];
        pcs=[pcs phases(i,2)*ones(1,phases(i,5)-phases(i,4)+1)];
    end
else
    res = permute(xce,[3 2 1]);
    q=squeeze(sum(sum(teste.^2,3),2));
    pcs=zeros(1,s(1));
end

[h95,r95,h,r]=ploton(t2,q,res,s(3),['k.-'],pcs,opt,alph,alpr,alph95,alpr95,s_sBn,axes1,axes2);

if phases(i,3) == s(1)-1
    plotoff(res,T2v,ssqres,T2vpred,ssqrespred,posTest,s(3),pcs(1),opt,alpoh,alpor,alpoh95,alpor95,axes3,axes4)
end

ph95=sum(h95)/length(h95);
ph=sum(h)/length(h);
pr95=sum(r95)/length(r95);
pr=sum(r)/length(r);
detect = r+h>0;
pdet=sum(detect)/length(detect);

