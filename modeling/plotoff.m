function [lima,limb,limar,limbr]=plotoff(resmod,hotelling,residuals,htest,rtest,pos,lotes,pc,opt,alpoh,alpor,alpoh95,alpor95,axes1,axes2)

% Plots the overall D-statistic and SPE values of the calibration batches using 
%   the leave-one-out cross-validated and theorical control limits. 
%
% [lima,limb,limar,limbr]=plotoff(resmod,hotelling,residuals,htest,rtest,pos,lotes,pc)                                             % call with standard parameters
%
% [lima,limb,limar,limbr]=plotoff(resmod,hotelling,residuals,htest,rtest,pos,lotes,pc,opt,alpoh,alpor,alpoh95,alpor95,axes1,axes2) % complete call
%
%
% INPUTS:
%
% resmod: (IxJxK) residuals in the calibration data, K(sampling times) 
%       x J(variables) x I(batches)
%
% hotelling: (Ix1) D-statistic of the calibration batches obtained in a
%    leave-one-out cross-validation, I(batches).
%
% residuals: (Ix1) SPE of the calibration batches obtained in a
%    leave-one-out cross-validation, I(batches).
%
% htest: (1x1) overall D-statistic for the test batch.
%
% rtest: (1x1) overall SPE for the test batch.
%
% pos: (1x1) position on the graph of the values estimated for the test batch.
%
% lotes: (1x1) number of calibration batches.
%
% pc: (1x1) number of principal components for the D-statistic.
%
% opt: boolean (1x1) 
%       true: plot results.
%       false: do not plot results.
%
% alpoh: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the D-statistic. 
%
% alpor: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the SPE.
%
% alpoh95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the D-statistic. 
%
% alpor95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the SPE.
%
% axes1: handle to the axes where the D-statistic chart is plotted.
%
% axes2: handle to the axes where the SPE chart is plotted.
%
%
% OUTPUTS:
%
% lima: D-statistic control limit at the suggested imposed significance level (alpha) for the 95% confidence 
%   level. 
%
% limb: D-statistic control limit at the suggested imposed significance level (alpha) for the 99% confidence 
%   level.
%
% limar: SPE control limit at the suggested imposed significance level (alpha) for the 95% confidence 
%   level. 
%
% limbr: SPE control limit at the suggested imposed significance level (alpha) for the 99% confidence 
%   level.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)
%            Jose Camacho Paez (josecamacho@ugr.es)
%           
% last modification: Oct/13
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
if nargin < 8, error('Incorrect number of input parameters. Please, check the help for further details.'); end;
sh = size(hotelling);
sr = size(residuals);
if sh ~= sr, error('The dimensions of Hotelling-T2 and SPE vectors do not coincide.'); end
if nargin < 9, opt = 1; end;

if nargin < 10, alpoh = 0.01; end;
if nargin < 11, alpor = 0.01; end;
if nargin < 12, alpoh95 = 0.05; end;
if nargin < 13, alpor95 = 0.05; end;

if nargin < 14 && opt, 
    figure;
    axes1 = axes; 
end;
if nargin < 15 && opt, 
    figure;
    axes2 = axes; 
end;

% D-statistic

s=size(resmod);
lclu = s(1);
if pos > s(1), lclu = pos; end

    if pc(1)~=0,
        lima =(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-alpoh95,pc,lotes-pc); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
        limb =(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-alpoh,pc,lotes-pc); % limite al 99% de conf.

        if opt,
            hot = hotelling./(limb.*ones(sh(1),1));
            axes(axes1)
            hold off
            cla
            bar(hot,'b');
            hold on
            bar(pos,htest/limb,'k');          
            plot(1:lclu,repmat(lima/limb,lclu,1),'r--');
            plot(1:lclu,repmat(limb/limb,lclu,1),'r');
            xlabel('Batches','FontSize', 12,'FontWeight','bold');
            ylabel('D-statistic','FontSize', 12,'FontWeight','bold');
            axis([0 lclu+1 0 max(1, max(max(hot),htest/limb))*1.05]);
        end
    else
        if opt,
            axes(axes1)
            hold off  
            plot(0,0)
        end
    end

% SPE 
    
   % Estimation of the SPE control limits following Jackson & Mudholkar's
   % approach
   %limar= spe_lim(resmod,alpor95); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
   %limbr= spe_lim(resmod,alpor); % limite al 99% de conf.  
    
   % Estimation of the SPE control limits following Box's approximation
   % approach
     E = unfold(permute(resmod,[3,2,1]),Inf);
     m = mean(sum(E.^2,2));
     v = var(sum(E.^2,2));
     limar = (v/(2*m))*chi2inv(1-alpor95,(2*m^2)/v);
     limbr = (v/(2*m))*chi2inv(1-alpor,(2*m^2)/v);
    
if opt,
    res = residuals./(limbr.*ones(sr(1),1));
    axes(axes2)
    hold off
    cla
    bar(res,'b');
    hold on;
    bar(pos,rtest/limbr,'k');
    if sr(1)>s(1), bar(sr(1),residuals(sr(1)),'k'); end
    plot(1:lclu,repmat(limar/limbr,lclu,1),'r--');
    plot(1:lclu,repmat(limbr/limbr,lclu,1),'r');
    xlabel('Batches','FontSize', 12,'FontWeight','bold');
    ylabel('Q-statistic','FontSize', 12,'FontWeight','bold');
    axis([0 lclu+1 0 max(1,max(max(res),rtest/limbr))*1.05]);
end