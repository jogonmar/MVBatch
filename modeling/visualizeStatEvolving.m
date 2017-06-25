function visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,label,plot_handles)
% 
% Visualizes the control charts of a multivariate statistic in its online version with its
% corresponding 95% and 99% control limits.
%
% visualizeStatEvolving(stat,lim95,lim99)                                                 % Minimum call to display the control chart for a multivariate statistic and the 95% and 99% control limits.
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv)                                 % Call to display the control chart for a multivariate statistic and the 95% and 99% control limits computed with theoretical distributions and adjusted by cross-validation).
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag)                     % As the previous call, with the option of normalizing the statistic using the control limit.
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,label)               % As previous call, but indicating the label to display on the y-axis
%
% visualizeStatEvolving(stat,lim95,lim99,lim95cv,lim99cv,scalingFlag,label,plot_handles)  % As previous call, but indicating handles to the figure where to display
%
% multiphaseFit(xini,phases,prep,flagCont)    % complete call 
%
% INPUTS
%
% stat: (KxI) online D or Q statistic values for I batches
%
% lim95: (Kx1) theoretical 95% control limit of the online D statistic 
%
% lim99: (Kx1) theoretical 99% control limit of the online D statistic 
%
% lim95cv: (Kx1) cross-validated 95% control limit of the online D
% statistic (if theoretical and cross-validated control limits need to be
% displayed)
%
% lim99cv: (Kx1) cross-validated 99% control limit of the online D statistic
% (if theoretical and cross-validated control limits need to be
% displayed)
%
% scalingFlag: (1x1) boolean to indicate whether the multivariate
% statistics must be normalized to the 99% control limits.
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

%% Paramter checking
assert(nargin >= 3, 'Incorrect number of input parameters');
if nargin < 4, lim95cv = []; end
if nargin < 5, lim99cv = []; end
if nargin < 6, scalingFlag = 0; end
if isempty(scalingFlag), scalingFlag = 0; end
assert(scalingFlag>=0 && scalingFlag<=1,'Wrong scaling parameter input. Note that this a boolean parameter to indicate whether multivariate statistics must be standardized to the control limits (1) or not (0)')
if nargin < 7, label = 'multivariate statistic'; end
if nargin < 8 
    figure;
    plot_handles = axes; 
end
    
s = size(stat);

if scalingFlag && isempty(lim99cv)
   stat = stat./repmat(lim99,1,s(2));
   lim95 = lim95./lim99;
   lim99 = lim99./lim99;
end

%% Main code
symbol = 'x';
if s(2) == 1, symbol = 'k.-'; end

cla(plot_handles);
axes(plot_handles);
hold off
plot(stat,symbol);
hold on

plot(real(lim95),'r--','Linewidth',1.2);
plot(real(lim99),'r','Linewidth',1.2);
if ~isempty(lim95cv) && ~isempty(lim99cv)
    plot(real(lim95cv),'g--');
    plot(real(lim99cv),'g');
end
v=axis;
axis([1 s(1) v(3) 2*max(lim99)]); 
xlabel('Sampling time','FontSize', 12,'FontWeight','bold');
ylabel(label,'FontSize', 12,'FontWeight','bold');