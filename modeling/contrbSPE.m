function [ovContrib,jContrib,kContrib] = contrbSPE(e,J,phases, axes1, axes2, axes3)

% Contribution plot to SPE.  
%
% y = ovContribSPE(e,J)       % Contribution plot to SPE without ploting phases
% y = ovContribSPE(e,J,phases)     % complete call
%
% INPUTS:
%
% e: (KJx1) array containing the residuals from a specific batch
%
% J: (1x1) number of process variables 
%
% phases: (nx1) number of phases/stages observed in the batch process (empty by default).
%
%
% OUTPUTS:
%
% ovovContrib: (JKx1) overall contribution to SPE
%
% jContrib: (Jx1) contribution to SPE per variable
%
% kContrib: (Kx1) contribution to SPE per batch time
%
% axes1: handle to the axes where the overall contribution to Hotelling-T2 is plotted.
%
% axes2: handle to the axes where the overall contribution to Hotelling-T2 per variable is plotted.
%
% axes3: handle to the axes where the overall contribution to Hotelling-T2 per batch time is plotted.
%
% codified by: Jose Gonzalez Martinez
% version: 1.0


% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, phases=[];end
if nargin < 4, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 5, 
    h2 = figure;
    axes2 = axes; 
end;
if nargin < 6, 
    h3 = figure;
    axes3 = axes; 
end;

% Initialization
K = size(e,1)/J;
jContrib = zeros(J,1);
kContrib = zeros(K,1);

% Calculation
ovContrib = (e.^2);

for j=1:J
   jContrib(j) = sum(ovContrib((j-1)*K+1:K*j),1);
end
for k=1:K
   kContrib(k)= sum(ovContrib(k:K:end,1));
end

   
axes(axes1);
bar(ovContrib); hold on;
for j=1:J
    line([j*K j*K],[min(ovContrib) max(ovContrib)],'LineStyle','--','Color','k');
end
xlabel('Variables','FontSize',16);
ylabel('Contribution to SPE ','FontSize',16);
set(h,'Color','white');
axis tight

axes(axes2);
bar(jContrib); set(h,'Color','white');axis tight;
xlabel('Original Process Variables','FontSize',16);
ylabel('Contribution to SPE','FontSize',16);

axes(axes3);
bar(kContrib); set(h,'Color','white');axis tight;
if not(isempty(phases))
   for p=1:numel(phases)
       line([phases(p) phases(p)],[min(kContrib) max(kContrib)],'LineStyle','--','Color','black');
   end
end
xlabel('Batch Time','FontSize',16);
ylabel('Contribution to SPE','FontSize',16);
