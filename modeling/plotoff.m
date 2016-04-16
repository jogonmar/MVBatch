function []=plotoff(resmod,hotelling,residuals,htest,rtest,pos,lotes,pc,opt,alpoh,alpor,alpoh95,alpor95,axes1,axes2)

% Plots the overall D-statistic and SPE values of the calibration batches using 
%   the leave-one-out cross-validated and theorical control limits. 
%
% [alph,alpr,alph95,alpr95]=plot_onOverallstats(resmod,hotelling,residuals,lotes,pc,opt) % call with standard parameters
%
% [alph,alpr,alph95,alpr95]=plot_onOverallstats(resmod,hotelling,residuals,lotes,tg
%   ,pc,opt,axes1,axes2) % complete call
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
% htest: (1x1) D-statistic of the test batch.
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
% alph: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the D-statistic. 
%
% alpr: suggested imposed significance level (alpha) for the 99% confidence 
%   limit in the SPE.
%
% alph95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the D-statistic. 
%
% alpr95: suggested imposed significance level (alpha) for the 95% confidence 
%   limit in the SPE.
%
%
% codified by: José M. Gonzalez-Martinez.
% version: 0.0
% last modification: 20/Oct/13.

% Parameters checking

if nargin < 7, error('The number of argument is not correct.'); end;
sh = size(hotelling);
sr = size(residuals);
if sh ~= sr, error('The dimensions of Hotelling-T2 and SPE vectors do not coincide.'); end
if nargin < 9, opt = 1; end;

if nargin < 10, alpoh = 0.01; end;
if nargin < 11, alpor = 0.01; end;
if nargin < 12, alpoh95 = 0.05; end;
if nargin < 13, alpor95 = 0.05; end;

if nargin < 14, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 15, 
    h2 = figure;
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
    ylabel('SPE','FontSize', 12,'FontWeight','bold');
    axis([0 lclu+1 0 max(1,max(max(res),rtest/limbr))*1.05]);
end