function [] = plot_onwarp(warpnoc,band,warptest,tg,nsamplesToPlot,plotallFlag,axes1)

% Design of the NOC Warping Information Control Chart. The original paper is Chemometrics and Intelligent
% Laboratory Systems, 127:210-217, 2013.
%
%  plot_onwarp(warpnoc,band)                                                    % call to build the NOCWI control chart.
%  plot_onwarp(warpnoc,band,warptest)                                           % call to build the NOCWI control chart and display the WI profile of a test batch with the default values.
%  plot_onwarp(warpnoc,band,warptest,tg,nsamplesToPlot,plotallFlag,axes1)       % complete call
%
% INPUTS:
%
% warpnoc: (KxI) calibration NOC warping profiles, K(time points) X J (batches)
%
% band: [Kx2] upper and lower boundaries constraining the search space in batch synchronization, K(time points)
%
% warptest: [Kx1] warping profile of the test batch, K(time points)
%
% tg: (1x1) marker to display on the control chart ('x' by default)
%
% nsamplesToPlot: (1x1) number of samples to display on the control chart (K by default)
%
% plotallFlag: (1x1) boolean to indicate whether the NOC WI profiles must be displayed (1 by default).
%
% axes1: handle to the axes where the warping profiles are plotted.
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: Aug/14
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

% Parameter checking
if nargin < 2, errodlg('Incorrect number of input parameters. Please, check the help for further details.');end
if nargin < 3, warptest = []; end   
if nargin < 4 || isempty(tg), tg = 'x';end
if nargin < 5 || isempty(nsamplesToPlot), nsamplesToPlot = size(warpnoc,1);end
if nargin < 6,  plotallFlag = 1;end
if nargin < 7, 
    figure;
    axes1 = axes; 
end;

s = size(warpnoc);
samplepoints = [1:s(1)]';
warpnocnew = warpnoc - repmat(samplepoints,1,s(2));
newband = zeros(s(1),2);

for i=1:s(1)
    newband(i,1) = band(i,1) - samplepoints(i);
    newband(i,2) = band(i,2) - samplepoints(i);
end
        
axes(axes1)
if  plotallFlag, plot(warpnocnew,tg,'MarkerSize',5); hold on; end
if ~isempty(warptest), plot(warptest(1:nsamplesToPlot) - repmat([1:nsamplesToPlot]',1,1),'kx-','MarkerSize',5); hold on;end

plot(1:s(1),newband,'r-'); hold on;
plot(1:s(1),ones(s(1),1),'k-');
v = axis;
v3=min(newband(:,2));
axis([1,s(1),min(v(3),v3)*1.15,max(v(4),v3)*1.15]);
xlabel('Sampling time (Reference)','FontSize', 12,'FontWeight','Bold');
ylabel('Sampling time (Test)','FontSize', 12,'FontWeight','Bold');
hold off;
    