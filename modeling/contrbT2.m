function [ovContrib,jContrib,kContrib] = contrbT2(x,T,P,J,phases,axes1,axes2,axes3)

% Contribution plot to Hotelling T2.  
%
% y = contrbT2(e,J)            % Contribution plot to Hotelling-T2 without ploting phases
% y = contrbT2(e,J,phases)     % complete call
%
% INPUTS:
%
% x: (KJx1) array containing the raw batch data for J process variables collected at K sampling points
%
% T: (IxA)  score matrix from the PCA model fitted on the training data set X.
%
% P: (KJxA) loading matrix from the PCA model fitted on the training data set X.
%
% J: (1x1) number of process variables 
%
% phases: (nx1) number of phases/stages observed in the batch process (empty by default).
%
% axes1: handle to the axes where the overall contribution to Hotelling-T2 is plotted.
%
% axes2: handle to the axes where the overall contribution to Hotelling-T2 per variable is plotted.
%
% axes3: handle to the axes where the overall contribution to Hotelling-T2 per batch time is plotted.
%
% OUTPUTS:
%
% ovovContrib: (JKx1) overall contribution to Hotelling-T2
%
% jContrib: (Jx1) contribution to Hotelling-T2 per variable
%
% kContrib: (Kx1) contribution to Hotelling-T2 per batch time
%
%
% codified by: Jose Gonzalez Martinez
% version: 1.0


% Parameters checking

if nargin < 4, error('Error in the number of arguments.'); end;
if nargin < 5, phases=[];end
if size(T,2) ~= size(P,2), error('Inconsistency in the number of latent variables'); end
if size(x,1) ~= size(P,1), error('Inconsistency in the number of process variables'); end
if nargin < 6, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 7, 
    h2 = figure;
    axes2 = axes; 
end;
if nargin < 8, 
    h3 = figure;
    axes3 = axes; 
end;

% Initialization
K = size(x,1)/J;
jContrib = zeros(J,1);
kContrib = zeros(K,1);

% Calculation
Sinv = inv(cov(T));  
ovContrib = real(x'*P*sqrt(Sinv)*P')';

for j=1:J
    jContrib(j,1) = sum(ovContrib((j-1)*K+1:K*j).^2);
end

for k=1:K
    kContrib(k,1)= sum(ovContrib(k:K:size(ovContrib,1)).^2);
end

% Contribution Plots

axes(axes1)
bar(ovContrib); hold on;
for j=1:J
    line([j*K j*K],[min(ovContrib) max(ovContrib)],'LineStyle','--','Color','k');
end
xlabel('Variables','FontSize',16);
ylabel('Contribution to Hotelling-T^2','FontSize',16);
set(h,'Color','white');
axis tight

axes(axes2)
bar(jContrib); set(h,'Color','white');axis tight;
xlabel('Original Process Variables','FontSize',16);
ylabel('Contribution to Hotelling-T^2','FontSize',16);

axes(axes3);
bar(kContrib); set(h,'Color','white');axis tight;
if not(isempty(phases))
   for p=1:numel(phases)
       line([phases(p) phases(p)],[min(kContrib) max(kContrib)],'LineStyle','--','Color','black');
   end
end
xlabel('Batch Time','FontSize',16);
ylabel('Contribution to Hotelling-T^2','FontSize',16);