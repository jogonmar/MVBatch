function [xrec] = reconstructX(X,T,P,pcs,meanv,stdv) 

% Batch trajectory forecast from the current time point k to the end of the
% batch using the Trimmed Score Regression. J. Chemometrics 2005; 19: 439–447
%
% [xrec] = reconstructX(X,T,P,pcs,meanv,stdv)        % complete call
%
% INPUTS:
%
% X: (kxJ) Two-way batch data matrix, k(observations) X J (variables)
%
% T: (IxH) scores where H stands for the maximum number of PCs.
%
% P: (JKxH) loadings where H stands for the maximum number of PCs.
%
% pc: (1x1) number of principal components.
%
% meanv: (KxJ) estimated average values 
%
% stdv: (KxJ) estimated standard deviation values 
%
% OUTPUTS:
%
% xrec: (1x1) TSR reconstruction of X with pc components.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: Apr/15
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

% Parameters checking
routine=dbstack;
assert (nargin >= 6, 'Error in the number of input parameters. Type ''help %s'' for more info.', routine(1).name);

% Validate dimensions of input data
Ht = size(T, 2);
Hp = size(P, 2);
A = size(T, 2);

assert (isequal(size(T,2), Hp), 'Dimension Error: number of PCs in the score matrix must be the same as in the loading matrix. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(P,2), Ht), 'Dimension Error: number of PCs in the loading matrix must be the same as in the score matrix. Type ''help %s'' for more info.', routine(1).name);
assert (A>0 && isequal(A,Hp) || isequal(A,Ht), 'The selected number of PCs must be lower or equal than the number of PCs extracted in the calibration model. Type ''help %s'' for more info.', routine(1).name);

% Initialization
l = size(X,1);
%Preprocess
xtest = (X - meanv(1:l,:))./stdv(1:l,:);
%Unfolding
xtestu = unfold(xtest,l-1);
jindb = 1:size(X,2)*l;
%TSR estimator
theta = cov(T);
theta_A = cov(T(:,1:pcs));
tauTSR = theta_A*P(jindb,1:pcs)'*P(jindb,1:pcs)*inv(P(jindb,1:pcs)'*P(jindb,:)*theta*P(jindb,:)'*P(jindb,1:pcs))*P(jindb,1:pcs)'*xtestu';
xpred = tauTSR'*P(:,1:pcs)';
% Reconstruction of the batch trajectory
xrec = fold([xtestu xpred(size(xtestu,2)+1:end)],1,size(meanv,1)-1);
xrec = xrec.*stdv+meanv;
