function [contribution] = contribution2Scores(X,T,P,pcs,indTo,indFrom)

% Contributions to scores between two groups of data. If you want to
% estimate the contribution plot of one group of samples to the average,
% indFrom must be set Inf.
%
% X: (KxJ) preprocess data matrix used for modelling
%
% T: (IxA) scores matrix
%
% P: (JxA) loading matrix
%
% pcs: (2x1) PCs taken into account for estimation
%
% indFrom: (nx1) indeces of the first group
%
% indTo: (nx1) indeces of the second group

% Checking

if nargin < 5, error('The number of parameters is incorrect');end
s = size(X);
if size(T,1) ~= s(1), error('The number of observations of the data and score matrices mismatches.'); end
if size(P,1) ~= s(2), error('The number of variables of the data and loading matrices mismatches.'); end
if size(pcs,1) ~= 2, error('Number of PCs unexpected');end
if indTo == Inf, indTo = 1:size(T,1);end
if nargin < 6, indFrom = 1:size(T,1); end
if indFrom == Inf, indFrom = 1:size(T,1);end

cov_inv=inv(cov(T));

Xto = mean(X(indTo,:),1);
Xfrom = mean(X(indFrom,:),1);

Tto = mean(T(indTo,pcs),1);
Tfrom = mean(T(indFrom,pcs),1);

contribution = zeros(4,1);

for j=1:size(P,1)
dist=0;
    for a=1:length(pcs)
        dist =  dist  + (((Tto(pcs(a))-Tfrom(pcs(a)))/sqrt(cov_inv(pcs(a),pcs(a))))^2)*P(j,pcs(a))^2;
    end 
    contribution(j) = (Xto(:,j)-Xfrom(:,j)).*sqrt(dist);
end


