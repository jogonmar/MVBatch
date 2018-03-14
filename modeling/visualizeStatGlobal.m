function visualizeStatGlobal(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,batchToHighlight,label,plot_handles)

% 
% Visualizes the control charts of a multivariate statistic in its offline version with its
% corresponding 95% and 99% control limits.
%
% visualizeStatEvolving(stat,lim95,lim99)                                                                      % Minimum call to display the control chart for a multivariate statistic and the 95% and 99% control limits.
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv)                                                      % Call to display the control chart for a multivariate statistic and the 95% and 99% control limits computed with theoretical distributions and adjusted by cross-validation).
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag)                                          % As the previous call, with the option of normalizing the statistic using the control limit.
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,batchToHighlight,label)                   % As previous call, but with the possibility to highlight one of values corresponding to a batch
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,batchToHighlight,label)                   % As previous call, but indicating the label to display on the y-axis
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,batchToHighlight,label,plot_handles)      % As previous call, but indicating handles to the figure where to display
%
% multiphaseFit(xini,phases,prep,flagCont)    % complete call 
%
% INPUTS
%
% stat: (Ix1) D or Q statistic values for I batches
%
% lim95: (1x1) theoretical 95% control limit of the online D statistic 
%
% lim99: (1x1) theoretical 99% control limit of the online D statistic 
%
% lim95cv: (1x1) cross-validated 95% control limit of the online D
% statistic (if theoretical and cross-validated control limits need to be
% displayed)
%
% lim99cv: (1x1) cross-validated 99% control limit of the online D statistic
% (if theoretical and cross-validated control limits need to be
% displayed)
%
% scalingFlag: (1x1) boolean to indicate whether the multivariate
% statistics must be normalized to the 99% control limits.
%
% batchToHighlight: (1x1) index (batch) whose statistic must be highlighted
% on the plot.
%
% label: string to display on y-axis
%
% plot_handles: (1x1) handles to the figure to display the control chart
%
% coded by: José M. González Martínez (jogonmar@gmail.com)
%
% Copyright (C) 2017  José M. González Martínez
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

%% Parameter checking
assert(nargin >=3, 'Incorrect number of input parameters');
s = size(stat);
if nargin < 4, lim95cv = []; end
if nargin < 5, lim99cv = []; end
if nargin < 6, scalingFlag = 0; end
if isempty(scalingFlag), scalingFlag = 0; end
assert(scalingFlag>=0 && scalingFlag<=1,'Wrong scaling parameter input. Note that this a boolean parameter to indicate whether multivariate statistics must be standardized to the control limits (1) or not (0)');
if nargin < 7, batchToHighlight = Inf; end
if isempty(batchToHighlight), batchToHighlight = Inf; end
if ~isinf(batchToHighlight) && (batchToHighlight < 1 || batchToHighlight > s(1)), error('Incorrect index to the batch to highylight. Ensure that you select one of the batches available.'); end
if nargin < 8, label = 'multivariate statistic'; end
if nargin < 9 
    figure;
    plot_handles = axes; 
end

if scalingFlag && isempty(lim99cv)
   stat = stat./repmat(lim99,1,s(2));
   lim95 = lim95./lim99;
   lim99 = lim99./lim99;
end

%% Main code
cla(plot_handles);
axes(plot_handles);
hold off
bar(stat,'b');
hold on;
if ~isinf(batchToHighlight), bar(batchToHighlight,stat(batchToHighlight),'k'); end

plot(0:s(1)+1,repmat(real(lim95),s(1)+2,1),'r--');
plot(0:s(1)+1,repmat(real(lim99),s(1)+2,1),'r');
if ~isempty(lim95cv) && ~isempty(lim99cv)
    plot(0:s(1)+1,repmat(real(lim95cv),s(1)+2,1),'g--');
    plot(0:s(1)+1,repmat(real(lim99cv),s(1)+2,1),'g');
end
v=axis;
axis([0 s(1)+1 v(3) max(v(4),1.2)]);
xlabel('Batches','FontSize', 12,'FontWeight','bold');
ylabel(label,'FontSize', 12,'FontWeight','bold');