function [P,T,mn,stnd,res,limdo95,limdo99,limqo95,limqo99,limd95,limd99,limq95,limq99] = multiphaseFit(xini,phases,prep,flagCont)

% Computes the multiphase model parameters for a calibration data set, the D and Q statistical values and its control limits at 95% and 99% confidence level. 
%
% multiphaseFit(xini,phases,prep)             % call without estimating contribution to statistics (by default)
%
% multiphaseFit(xini,phases,prep,flagCont)    % complete call 
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
% flagCont: (1x1) boolean to indicate whether the contributions to D and Q
% statistics must be computed.
%
% OUTPUTS 
%
% P: {phases x 1} cell array of [J_phase x A_phase] matrices of loadings for all phases.
%
% T: {phases x 1}cell array of [I_phase x A_phase] matrices of scores for all phases. 
%
% mn: [KxJ] matrix of averages.
%
% stnd: [KxJ] matrix of standard deviations.
%
% res: (IxJxK) residuals from the calibration data set, K(sampling times) 
%       x J(variables) x I(batches)
%
% limdo95: (Kx1) theoretical 95% control limit of the online D statistic 
%
% limdo99: (Kx1) theoretical 99% control limit of the online D statistic 
%
% limqo95: (Kx1) theoretical 95% control limit of the online Q statistic 
%
% limqo99: (Kx1) theoretical 99% control limit of the online Q statistic  
%
% limd95: (1x1) theoreticald 95% control limit of the offline D statistic (only computed for batch-wise modeling) 
%
% limd99: (1x1) theoretical 99% control limit of the offline D statistic (only computed for batch-wise modeling)
%
% limq95: (1x1) theoretical 95% control limit of the offline Q statistic (only computed for batch-wise modeling)
%
% limq99: (1x1) theoretical 99% control limit of the offline Q statistic (only computed for batch-wise modeling)
%
% coded by: José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2017  José M. González Martínez
% Copyright (C) 2017  Jose Camacho Paez, University of Granada, Granada
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

%% Parameters checking
routine=dbstack;
assert (nargin >=2, 'Error in the number of input parameters. Type ''help %s'' for more info.', routine(1).name);
assert (ndims(xini)==3, 'Incorrect number of dimensions of xini.');
s = size(xini);
assert (isempty(find(s<1)), 'Incorrect number of dimensions of xini.');
if nargin <= 3, prep = 2; end;
assert (prep>=0 && prep<=5, 'Incorrect content of prep.');
if nargin < 4, flagCont = 0; end
assert(flagCont == 0 || flagCont == 1,'Incorrect flag value to compute contribution. It must be either 0 (no compute) or 1 (compute).'); 


% Initialization
sp     = size(phases);
res    = zeros(s(3),s(2),s(1));
T      = cell(sp(1),1);
P      = cell(sp(1),1);
pcs    = zeros(s(1),1);
alpdo95 = 0.05;
alpdo99 = 0.01;
alpqo95 = 0.05;
alpqo99 = 0.01;
limdo95 = nan(s(1),1); 
limqo95 = nan(s(1),1);
limdo99 = nan(s(1),1);
limqo99 = nan(s(1),1);
alpd95 = 0.05;
alpd99 = 0.01;
alpq95 = 0.05;
alpq99 = 0.01;
limd95 = NaN; 
limq95 = NaN;
limd99 = NaN;
limq99 = NaN;

%% Main code
[xce,mn,stnd] = preprocess3D(xini,prep);

for ph=1:sp(1),
    ind=(max(phases(ph,4)-phases(ph,3),1):phases(ph,5));
    ind_ini=find(ind==phases(ph,4));

    xu=unfold(xce(ind,:,:),phases(ph,3));
    [U,S,V] = svd(xu,'econ');
    T{ph,1} = U*S;
    P{ph,1} = V;        
    p = P{ph,1}(:,1:phases(ph,2));
    t = T{ph,1}(:,1:phases(ph,2));
    resa=xu-t*p';
    resb=fold(resa,s(3),phases(ph,3));
    res(:,:,ind(ind_ini:end))=permute(resb(ind_ini:end,:,:),[3 2 1]);

    if ind_ini<phases(ph,3)+1,
    menor_en=phases(ph,3)-ind_ini;
    theta = cov(T{ph,1});
    theta_A = cov(t);
        for j=0:menor_en,
            indb=1:s(2)*(ind_ini+j);
            ind2=phases(ph,4)+j;
            % IMPUTATION USING TSR
            t_t = theta_A*p(indb,:)'*p(indb,:)*inv(p(indb,:)'*P{ph,1}(indb,:)*theta*P{ph,1}(indb,:)'*p(indb,:))*p(indb,:)'*xu(1:s(3),indb)';
            res(:,:,ind2)=permute(xce(ind2,:,:),[3 2 1])-t_t'*p(indb(end-s(2)+1:end),:)';   
        end        
    end

    % Estimate the control limits
    pcs(phases(ph,4):phases(ph,5)) = phases(ph,2)*ones(1,phases(ph,5)-phases(ph,4)+1); 
end
if flagCont
    % Compute the theoretical control limits for on-line monitoring
    for i=1:s(1)            
        limdo95(i) = (pcs(i)*(s(3)*s(3)-1)/(s(3)*(s(3)-pcs(i))))*finv(1-alpdo95,pcs(i),s(3)-pcs(i)); 
        limdo99(i) = (pcs(i)*(s(3)*s(3)-1)/(s(3)*(s(3)-pcs(i))))*finv(1-alpdo99,pcs(i),s(3)-pcs(i)); 
        limqo95(i)= spe_lim_box(res(:,:,i),alpqo95); 
        limqo99(i)= spe_lim_box(res(:,:,i),alpqo99); 
    end
    % Compute the theoretical control limits for off-line monitoring if
    % model is batch-wise
    if sp(1)==1 && phases(1,3)== s(1)-1
        limd95 = (pcs(1)*(s(3)*s(3)-1)/(s(3)*(s(3)-pcs(1))))*finv(1-alpd95,pcs(1),s(3)-pcs(1)); 
        limd99 = (pcs(1)*(s(3)*s(3)-1)/(s(3)*(s(3)-pcs(1))))*finv(1-alpd99,pcs(1),s(3)-pcs(1)); 
        limq95 = spe_lim_box(resa,alpq95); 
        limq99 = spe_lim_box(resa,alpq99); 
    end
end