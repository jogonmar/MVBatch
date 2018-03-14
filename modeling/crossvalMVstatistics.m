function [cvevolD,cvevolQ,cvD,cvQ,limdo95cv,limdo99cv,limqo95cv,limqo99cv,alpdo95cv,alpdo99cv,alpqo95cv,alpqo99cv,limd95cv,limd99cv,limq95cv,limq99cv,alpd95cv,alpd99cv,alpq95cv,alpq99cv] = crossvalMVstatistics(xini,phases,prep)

% Computes D-statistic and SPE values for the calibration batches using leave-one-out cross-validation. 
%
% [alph,alpr,alph95,alpr95]=plot_distcv2(xini, phases, prep)            % complete call 
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
%
% OUTPUTS (class parameters):
%
% cvevolD: (Kx1) cross-validated online D-statistic values for the test batch.
%
% cvevolQ: (Kx1) cross-validated online Q-statistic values for the test batch.
%
% cvD: (Ix1) cross-validated overall D-statistic values for the test batch (if batch-wise modeling).
%
% cvQ: (Ix1) cross-validated overall Q-statistic values for the test batch (if batch-wise modeling).
%
% limdo95cv: (Kx1) cross-validated 95% control limit of the D statistic for online applications.
%
% limdo99cv: (Kx1) cross-validated 99% control limit of the D statistic for online applications.
%
% limqo95cv: (Kx1) cross-validated 95% control limit of the Q statistic for online applications.
%
% limqo99cv: (Kx1) cross-validated 99% control limit of the Q statistic ffor online applications.
%
% alpdo95cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the online D statistic. 
%
% alpdo99cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the online D statistic.
%
% alpqo95cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the online Q statistic. 
%
% alpqo99cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the online Q statistic.
%
% limd95cv: (1x1) cross-validated 95% control limit of the D statistic for offline applications.
%
% limd99cv: (1x1) cross-validated 99% control limit of the D statistic for offline applications.
%
% limq95cv: (1x1) cross-validated 95% control limit of the Q statistic for offline applications.
%
% limq99cv: (1x1) cross-validated 99% control limit of the Q statistic for offline applications. 
%
% alpd95cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the offline D statistic. 
%
% alpd99cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the offline D statistic.
%
% alpq95cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the offline Q statistic. 
%
% alpq99cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the offline Q statistic.
%
% 
% coded by:  José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2017 José M. González Martínez
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
assert (nargin >= 2, 'Error in the number of input parameters. Type ''help %s'' for more info.', routine(1).name);
assert (ndims(xini)==3, 'Incorrect number of dimensions of xini.');
s = size(xini);
assert (isempty(find(s<1)), 'Incorrect number of dimensions of xini.');
if nargin < 4, prep = 2; end;
assert (prep>=0 && prep<=5, 'Incorrect content of prep.');
sp = size(phases);
% Intialization
cvevolD = nan(s(1),s(3));
cvevolQ = nan(s(1),s(3));
cvD     = nan(s(3),1);
cvQ     = nan(s(3),1);
% Online control limits
limdo95cv = nan(s(1),1);
limdo99cv = nan(s(1),1);
limqo95cv = nan(s(1),1);
limqo99cv = nan(s(1),1);
% Offline control limits
limd95cv = nan(s(3),1);
limd99cv = nan(s(3),1);
limq95cv = nan(s(3),1);
limq99cv = nan(s(3),1);
% Imposed online levels
alpd95cv = 0.05;
alpd99cv = 0.05;
alpq95cv = 0.01;
alpq99cv = 0.01;
% Imposed offline levels
alpdo95cv = 0.05;
alpdo99cv = 0.05;
alpqo95cv = 0.01;
alpqo99cv = 0.01;

pcs     = [];


%% Main code

h = waitbar(0/s(3),'Leave-one batch-out cross-validation - Initializing','Name','Cross-validation of the monitoring system');
try
    [po,to,~,~,reson,limdo95,limdo99,limqo95,limqo99,limd95,limd99,limq95,limq99] = multiphaseFit(xini,phases,prep,1);

    for o=1:s(3)
        test=xini(:,:,o);
        xini2=xini(:,:,[1:o-1 o+1:s(3)]);

        % Calibrate model
        [p,t,av,stdn] = multiphaseFit(xini2,phases,prep);

        % Project test onto latent structure
        [evolD,evolQ,D,Q] = multiphaseProjection(xini2,test,phases,p,t,av,stdn);
        cvevolD(:,o) = evolD;
        cvevolQ(:,o) = evolQ;
        cvD(o)     = D;
        cvQ(o)     = Q;
        waitbar(o/s(3),h,sprintf('Leave-one batch-out cross-validation - Batch #%d',o),'Name','Cross-validation of the monitoring system');
    end 
catch err
   error(err.message); 
   close(h);
   return;
end
close(h);

% Estimate the control limits
for p=1:sp(1), pcs = [pcs; phases(p,2)*ones(1,phases(p,5)-phases(p,4)+1)]; end

% Cross-validate control limits
% On-line 
try
h = waitbar(1/2,'Computing control limits for instantaneous D- and Q-statistics','Name','Cross-validation of the monitoring system');
[limdo95cv,limdo99cv,limqo95cv,limqo99cv,alpdo95cv,alpdo99cv,alpqo95cv,alpqo99cv] = crossval_limits(reson,cvevolD,cvevolQ,limdo95,limdo99,limqo95,limqo99,s(3),pcs);
waitbar(1,h,'Computing online control limits for instantaneous D- and Q-statistics','Name','Cross-validation of the monitoring system');
catch err
   error(err.message); 
   close(h);
   return;
end
close(h);
% Off-line
try
    if sp(1) == 1 && phases(1,3) == s(1)-1
        xce = preprocess3D(xini,prep);
        xu=unfold(xce,phases(1,3));
        resoff = xu-to{1}(:,1:phases(1,2))*po{1}(:,1:phases(1,2))';
        h = waitbar(1/2,'Computing control limits for overall D- and Q-statistics','Name','Cross-validation of the monitoring system');
        [limd95cv,limd99cv,limq95cv,limq99cv,alpd95cv,alpd99cv,alpq95cv,alpq99cv] = crossval_limits(resoff,cvD,cvQ,limd95,limd99,limq95,limq99,s(3),repmat(pcs(1),s(3),1));
        waitbar(1,h,'Computing online control limits for overall D- and Q-statistics','Name','Cross-validation of the monitoring system');
    end
catch err
    error(err.message);
    close(h);
    return;
end
close(h);


end % CVMonitorParameters
