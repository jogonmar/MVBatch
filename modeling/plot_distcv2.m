function [alph,alpr,alph95,alpr95,alpoh,alpor,alpoh95,alpor95,overallOffMon]=plot_distcv2(xini, phases, prep, opt, axes1, axes2,axes3,axes4)

% Computes the D-statistic and SPE values of the calibration batches using 
%   leave-one-out cross-validation. 
%
% [alph,alpr,alph95,alpr95]=plot_distcv(xini, phases, prep, opt) 
%    % call with standard parameters
%
% [alph,alpr,alph95,alpr95,alpoh,alpor,alpoh95,alpor95,overallOffMon]=plot_distcv(xini, phases, prep, opt, axes1, 
%   axes2,axes3,axes4) % complete call
%
%
% INPUTS:
%
% xini: (KxJxI) three-way batch data matrix for calibration, K(sampling times) 
%       x J(variables) x I(batches)
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
% axes1: handle to the axes where the online D-statistic chart is plotted.
%
% axes2: handle to the axes where the online SPE chart is plotted.
%
% axes3: handle to the axes where the offline D-statistic chart is plotted.
%
% axes4: handle to the axes where the offline SPE chart is plotted.
%
%
% OUTPUTS:
%
% alph: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the D-statistic. 
%
% alpr: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the SPE.
%
% alph95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the D-statistic. 
%
% alpr95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the SPE.
%
% alpoh: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the overall D-statistic. 
%
% alpor: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the overall SPE.
%
% alpoh95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the overall D-statistic. 
%
% alpor95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the overall SPE.
%
% 
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: 1/Aug/14
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

if nargin < 2, error('Numero de argumentos erroneos.'); end;

if ndims(xini)~=3, error('Incorrect number of dimensions of xini.'); end;
s = size(xini);
if find(s<1), error('Incorrect content of xini.'); end;

if ndims(phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp=size(phases);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of phases.'); end;
if find(phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(phases(:,4:5)<1), error('Incorrect content of phases.'); end;
if find(phases(:,3:5)>s(1)), error('Incorrect content of phases.'); end;

if nargin < 3, prep = 2; end;
if nargin < 4, opt = 1; end;

if (prep<0||prep>5), error('Incorrect value of prep.'); end;

if nargin < 5, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 6, 
    h2 = figure;
    axes2 = axes; 
end;

if nargin < 7, 
    h3 = figure;
    axes3 = axes; 
end;
if nargin < 8, 
    h4 = figure;
    axes4 = axes; 
end;


% Main code

[xce,av,sta] = preprocess3D(xini,prep);
    
q=[];
res=zeros(s(3),s(2),s(1));
pcs=[];
for i=1:sp(1),
    ind=(max(phases(i,4)-phases(i,3),1):phases(i,5));
    ind_ini=find(ind==phases(i,4));
    xu=unfold(xce(ind,:,:),phases(i,3));
    [U,S,V] = svd(xu,'econ');
    tAll = U*S;
    pAll = V;
    p = pAll(:,1:phases(i,2));
    t = tAll(:,1:phases(i,2));
        
    resa=xu-t*p';
    resb=fold(resa,s(3),phases(i,3));
    res(:,:,ind(ind_ini:end))=permute(resb(ind_ini:end,:,:),[3 2 1]);

     if ind_ini<phases(i,3)+1,
        menor_en=phases(i,3)-ind_ini;
        for j=0:menor_en,
            indb=1:s(2)*(ind_ini+j);
            ind2=phases(i,4)+j;
            t_t = inv(p(indb,:)'*p(indb,:))*p(indb,:)'*xu(1:s(3),indb)';
            res(:,:,ind2)=permute(xce(ind2,:,:),[3 2 1])-t_t'*p(indb(end-s(2)+1:end),:)';   
            q=[q ;sum((permute(xce(ind2,:,:),[3 2 1])-t_t'*p(indb(end-s(2)+1:end),:)').^2,2)'];
        end        
    end
    
    ssqres = sum(resa.^2,2)';
    qb=[];
    for o=1:s(3):length(ssqres),
        qb=[qb;ssqres(o:o+s(3)-1)];
    end
    q=[q;qb];
    pcs=[pcs phases(i,2)*ones(1,phases(i,5)-phases(i,4)+1)];
    
    %% Code for post-batch off-line process monitoring

    if phases(i,3) == s(1)-1
        Sinv = inv(cov(t));
        T2v = zeros(s(3),1);
        for j=1:s(3)
            T2v(j) = (t(j,:))*Sinv*(t(j,:))';
        end
    end
end

tcv=[];
qcv=[];
for o=1:s(3),
    test=xini(:,:,o);
    xini2=xini(:,:,[1:o-1 o+1:s(3)]);
    tcvb=[];
    qcvb=[];
    for i=1:sp(1),
        ind=(max(phases(i,4)-phases(i,3),1):phases(i,5));
        ind_ini=find(ind==phases(i,4));
        
        [xce,av,sta] = preprocess3D(xini2(ind,:,:),prep);
        xu=unfold(xce,phases(i,3));
        [U,S,V] = svd(xu,'econ');
        tAll = U*S;
        pAll = V;
        p = pAll(:,1:phases(i,2));
        t = tAll(:,1:phases(i,2));
        teste = (test(ind,:,:)-av)./sta;
        testu=unfold(teste,phases(i,3));
        tpred = testu*p;
        resb=testu-tpred*p';
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
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*xu(1:(s(3)-1),jindb)';
                cov_inv=inv(cov(t_t'));
                t_t = theta_A*p(jindb,:)'*p(jindb,:)*inv(p(jindb,:)'*pAll(jindb,:)*theta*pAll(jindb,:)'*p(jindb,:))*p(jindb,:)'*testu(1,jindb)';
                % CALCULATE STATISTICS
                tcvb=[tcvb;t_t'*cov_inv*t_t];
                qcvb=[qcvb ;sum((permute(teste(jind2,:,:),[3 2 1])-t_t'*p(jindb(end-s(2)+1:end),:)').^2)];
            end        
        end
        ssc=size(t);
        j=1;
        for u=1:(s(3)-1):ssc(1),
            sc_model = t(u:u+s(3)-2,:);
            cov_inv = inv(cov(sc_model));
            tcvb=[tcvb;tpred(j,:)*cov_inv*tpred(j,:)'];
            j=j+1;
        end
    
        qcvb=[qcvb;squeeze(sum(sum(resb((phases(i,3)+1):end,:,:).^2,3),2))];
    end
    tcv=[tcv tcvb];
    qcv=[qcv qcvb];
end

alpoh95 = []; lima = [];limacv = [];
alpoh = []; limb  =[];limbcv  =[];
alpor95 = []; limar = []; limarcv = [];
alpor = []; limbr = []; limbrcv = [];
    

[alph,alpr,alph95,alpr95]=plotcv2(res,tcv,qcv,s(3),['x'],pcs,opt,axes1,axes2);
 
cla(axes3);
 cla(axes4);
set(axes3,'Visible','off');
 set(axes4,'Visible','off');
if phases(i,3) == s(1)-1
    [alpoh,alpor,alpoh95,alpor95]=plotOverall(resa,T2v,ssqres,s(3),pcs(1),opt,axes3,axes4);
     set(axes3,'Visible','on');
     set(axes4,'Visible','on');
end

overallOffMon = [];
% Save statistics from the overall off-line process monitoring in a
% structure
if phases(i,3) == s(1)-1
    overallOffMon.lima = lima;
    overallOffMon.limb = limb;
    overallOffMon.limar = limar;
    overallOffMon.limbr = limbr;
    overallOffMon.T2 = T2v;
    overallOffMon.SPE = ssqres';
end
