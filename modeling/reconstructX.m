function [xrec] = reconstructX(x,t,p,pcs,M,S) 

l = size(x,1);
 
%Preprocess
 xtest = (x - M(1:l,:))./S(1:l,:);
%Unfolding
 xtestu = unfold(xtest,l-1);
 jindb = 1:size(x,2)*l;

%TSR estimator
theta = cov(t);
theta_A = cov(t(:,1:pcs));

tauTSR = theta_A*p(jindb,1:pcs)'*p(jindb,1:pcs)*inv(p(jindb,1:pcs)'*p(jindb,:)*theta*p(jindb,:)'*p(jindb,1:pcs))*p(jindb,1:pcs)'*xtestu';
xpred = tauTSR'*p(:,1:pcs)';

% Reconstruction of the batch trajectory
s = size(M);
xrec2 = zeros(s);
xrec = fold([xtestu xpred(size(xtestu,2)+1:end)],1,size(M,1)-1);

xrec = xrec.*S+M;
