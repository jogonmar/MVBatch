function [limd95cv,limd99cv,limq95cv,limq99cv,alpd95cv,alpd99cv,alpq95cv,alpq99cv] = crossval_limits(resmod,cvD,cvQ,limd95,limd99,limq95,limq99,nbatches,pcs)

% Cross-validation of the control limimts for D and Q statistics using leave-one-out approach. 
%
% [limd95cv,limd99cv,limq95cv,limq99cv,alpd95cv,alpd99cv,alpq95cv,alpq99cv] = crossval_limits(cvD,cvQ,limd95,limd99,limq95,limq99,nbatches,pcs) % complete call
%
% INPUTS:
%
% resmod: [online version](IxJxK) residuals from the calibration data set, K(sampling times) x J(variables) x I(batches)
%         [offline version - batcvh-wise](IxKJ) residuals in the calibration data, K(sampling times) x J(variables) x I(batches)
%
% cvD:    [online version] (KxI) K online D-statistic values for I test batches.
%         [offline version - batcvh-wise] (Ix1) offline D-statistic values for I test batches.
%
% cvQ:    [online version] (KxI) K online Q-statistic values for I test batches.
%         [offline version - batcvh-wise] (Ix1) offline Q-statistic values of I test batches.
%
% limd95: [online version] (Kx1) 95% control limit of the D statistic 
%         [offline version - batcvh-wise] (1x1) 95% control limit of the D statistic
%
% limd99: [online version] (Kx1) 99% control limit of the D statistic 
%         [offline version - batcvh-wise] (1x1) 99% control limit of the D statistic
%
% limq95: [online version] (Kx1) 95% control limit of the Q statistic 
%         [offline version - batcvh-wise] (1x1) 95% control limit of the Q statistic
%
% limq99: [online version] (Kx1) 99% control limit of the Q statistic 
%         [offline version - batcvh-wise] (1x1) 99% control limit of the Q statistic
%
% nbatches:  (1x1) number of calibration batches.
%
% pc: [online version] (Kx1) number of principal components for the D-statistic over the batch run.
%     [offline version - batcvh-wise] (Ix1) number of principal components for the D-statistic across batches.
%
% OUTPUTS:
%
% limd95cv: (Kx1) or (1x1) cross-validated 95% control limit of the D statistic for online and offline applications, respectively.
%
% limd99cv: (Kx1) or (1x1) cross-validated 99% control limit of the D statistic for online and offline applications, respectively. 
%
% limq95cv: (Kx1) or (1x1) cross-validated 95% control limit of the Q statistic for online and offline applications, respectively.
%
% limq99cv: (Kx1) or (1x1) cross-validated 99% control limit of the Q statistic for online and offline applications, respectively. 
%
% alpd95cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the D statistic. 
%
% alpd99cv: (1x1) cross-validated significance level (alpha) for the 99% control limit of the D statistic.
%
% alpq95cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the Q statistic. 
%
% alpq99cv: (1x1) cross-validated significance level (alpha) for the 95% control limit of the Q statistic.
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

% Parameters checking
routine=dbstack;
assert (nargin == 9, 'Error in the number of input parameters. Type ''help %s'' for more info.', routine(1).name);

% Initialization
alpd95 = 0.05;
alpq95 = 0.05;
alpd99 = 0.01;
alpq99 = 0.01;
s=size(cvD);
limd95cv = nan(s(1),1);
limd99cv = nan(s(1),1);
limq95cv = nan(s(1),1);
limq99cv = nan(s(1),1);

if s(1) == size(cvD,1) && length(limd95)==1
   postbatch = 1; 
else
   postbatch = 0;
end

% Theoretical limits
if postbatch
    limd95 = repmat(limd95,s(1),1);
    limd99 = repmat(limd99,s(1),1);
end

%% D-statistic

% Experimental limit
cvevolD1way = reshape(cvD./repmat(limd95,1,size(cvD,2)),s(1)*s(2),1);
cvevolD1waySorted = sort(cvevolD1way);

lev = round(s(1)*s(2)*(1-alpd95));
ind = cvevolD1waySorted(lev);
ind2 = mod(find(cvevolD1way==ind,1),s(1));
if ~ind2, ind2 = s(1); end;
alpd95cv=1-fcdf(((nbatches*(nbatches-pcs(ind2)))/(pcs(ind2)*(nbatches*nbatches-1)))*(limd95(ind2) * ind),pcs(ind2),nbatches-pcs(ind2));

cvevolD1way = reshape(cvD./repmat(limd99,1,size(cvD,2)),s(1)*s(2),1);
cvevolD1waySorted = sort(cvevolD1way);

lev = round(s(1)*s(2)*(1-alpd99));
ind = cvevolD1waySorted(lev);
ind2 = mod(find(cvevolD1way==ind,1),s(1));
if ~ind2, ind2 = s(1); end;

alpd99cv=1-fcdf(((nbatches*(nbatches-pcs(ind2)))/(pcs(ind2)*(nbatches*nbatches-1)))*(limd99(ind2) * ind),pcs(ind2),nbatches-pcs(ind2));

if postbatch
    limd95cv = (pcs(1)*(nbatches*nbatches-1)/(nbatches*(nbatches-pcs(1))))*finv(1-alpd95cv,pcs(1),nbatches-pcs(1)); 
    limd99cv = (pcs(1)*(nbatches*nbatches-1)/(nbatches*(nbatches-pcs(1))))*finv(1-alpd99cv,pcs(1),nbatches-pcs(1)); 
else
    for i=1:s(1)
        limd95cv(i) = (pcs(i)*(nbatches*nbatches-1)/(nbatches*(nbatches-pcs(i))))*finv(1-alpd95cv,pcs(i),nbatches-pcs(i)); 
        limd99cv(i) = (pcs(i)*(nbatches*nbatches-1)/(nbatches*(nbatches-pcs(i))))*finv(1-alpd99cv,pcs(i),nbatches-pcs(i)); 
    end    
end

%% SPE 

% Identify the value close to the alpha % for the imposed level at 95%
cvevolQ95way = reshape(cvQ./repmat(limq95,1,size(cvQ,2)),s(1)*s(2),1);
cvevolQ95waySorted = sort(cvevolQ95way);

lev = round(s(1)*s(2)*(1-alpq95));
ind95 = cvevolQ95waySorted(lev);
ind95b = mod(find(cvevolQ95way==ind95,1),s(1));
if ~ind95b, ind95b = s(1); end;

% Identify the value close to the alpha % for the imposed level at 99%
cvevolQ99way = reshape(cvQ./repmat(limq99,1,size(cvQ,2)),s(1)*s(2),1);
cvevolQ99waySorted = sort(cvevolQ99way);

lev = round(s(1)*s(2)*(1-alpq99));
ind99 = cvevolQ99waySorted(lev);
ind99b = mod(find(cvevolQ99way==ind99,1),s(1));
if ~ind99b, ind99b = s(1); end

% Adjust of the confidence level to meet the imposed 95% and 99% confidence
% level following Box's approximation. If this fails, use Jackson & Mudholkar's approach
if postbatch
   alpq95cv = spe_pvalue_box(resmod,(limq95 * ind95)); limq95cv = spe_lim_box(resmod,alpq95cv);
   alpq99cv = spe_pvalue_box(resmod,(limq99 * ind99)); limq99cv = spe_lim_box(resmod,alpq99cv);
   if alpq95cv==0 || alpq99cv==0
       alpq95cv = spe_pvalue(resmod,(limq95 * ind95)); limq95cv = spe_lim(resmod,alpq95cv); 
       alpq99cv = spe_pvalue(resmod,(limq99 * ind99)); limq99cv = spe_lim(resmod,alpq99cv); 
   end
else
   alpq95cv =spe_pvalue_box(resmod(:,:,ind95b),(limq95(ind95b) * ind95));
   alpq99cv =spe_pvalue_box(resmod(:,:,ind99b),(limq99(ind99b) * ind99));
   
   for i=1:s(1)
         limq95cv(i)= spe_lim_box(resmod(:,:,i),alpq95cv); 
         limq99cv(i)= spe_lim_box(resmod(:,:,i),alpq99cv);  
   end
end  

end
