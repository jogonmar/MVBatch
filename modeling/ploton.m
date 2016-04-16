function [hi95,ri95,hi,ri]=ploton(hotelling,residuals,resmod,lotes,tg,pc,opt,alph,alpr,alph95,alpr95,s_sBn,axes1,axes2)

% Plots the D-statistic and SPE charts for on-line monitoring. 
%
% [hi95,ri95,hi,ri]=ploteaon(hotelling,residuals,resmod,lotes,tg,pc,opt) 
%   % call with standard parameters
%
% [hi95,ri95,hi,ri]=ploteaon(hotelling,residuals,resmod,lotes,tg,pc,opt,
%   alph,alpr,alph95,alpr95) % output in MATLAB console
%
% [hi95,ri95,hi,ri]=ploteaon(hotelling,residuals,resmod,lotes,tg,pc,opt,
%   alph,alpr,alph95,alpr95,axes1,axes2) % complete call
%
%
% INPUTS:
%
% hotelling: (Kx1) D-statistic of the test batch, K(sampling times).
%
% residuals: (Kx1) SPE of the test batch, K(sampling times).
%
% resmod: (IxJxK) residuals in the calibration data, K(sampling times) 
%       x J(variables) x I(batches)
%
% lotes: (1x1) number of calibration batches.
%
% tg: string with the plot options for the statistics of the test batch.
%
% pc: (1x1) number of principal components for the D-statistic.
%
% opt: boolean (1x1) 
%       true: plot results.
%       false: do not plot results.
%
% alph: imposed significance level (alpha) for the 99% confidence limit in the 
%   D-statistic. 
%
% alpr: imposed significance level (alpha) for the 99% confidence limit in the 
%   SPE.
%
% alph95: imposed significance level (alpha) for the 95% confidence limit in the 
%   D-statistic. 
%
% alpr95: imposed significance level (alpha) for the 95% confidence limit in the 
%   SPE.
%
% s_sBn: number of synchronized samples to be monitored (default: Inf, i.e.
% all).
%
% axes1: handle to the axes where the D-statistic chart is plotted.
%
% axes2: handle to the axes where the SPE chart is plotted.
%
%
% OUTPUTS:
%
% hi95: (1xK) vector with 1s in the sampling times where the D-statistic for 
%   the test batch exceeds the 95% control limit and 0s in the rest.
%
% ri95: (1xK) vector with 1s in the sampling times where the SPE for the test
%   batch exceeds the 95% control limit and 0s in the rest.
%
% hi: (1xK) vector with 1s in the sampling times where the D-statistic for 
%   the test batch exceeds the 99% control limit and 0s in the rest.
%
% ri: (1xK) vector with 1s in the sampling times where the SPE for the test
%   batch exceeds the 99% control limit and 0s in the rest.
%
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 13/Dic/11 by José M. González Martínez. The input
% parameter 's_sBn' has been added to control the number of samples to be
% plotted.

% Parameters checking

if nargin < 6, error('Numero de argumentos erroneos.'); end;
if nargin < 7, opt = 1; end;
if nargin < 8, alph = 0.01; end;
if nargin < 9, alpr = 0.01; end;
if nargin < 10, alph95 = 0.05; end;
if nargin < 11, alpr95 = 0.05; end;
if nargin < 12, s_sBn = Inf; end ;
if nargin < 13
    h = figure;
    axes1 = axes; 
end;
if nargin < 14, 
    h2 = figure;
    axes2 = axes; 
end;

if s_sBn == Inf, s_sBn = size(hotelling);end

% D-statistic
hi95 = [];
hi = [];

s=size(residuals);
out = false;

lclu = s(1);

if pc(1)~=0,
    lima=[];
    limb=[];
    for i=1:s(1),
        lima(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alph95,pc(i),lotes-pc(i)); 
        limb(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alph,pc(i),lotes-pc(i));
    end

    hot = hotelling./(limb'*ones(1,s(2)));
     if opt,
        axes(axes1)
        hold off
        plot(hot(1:s_sBn),tg)
        hold on;

        xlabel('Sampling time','FontSize', 12,'FontWeight','Bold');
        ylabel('D-statistic','FontSize', 12,'FontWeight','Bold');

        lim95=real(lima)./real(limb);
        plot(1:s(1),lim95,'r--');
        plot(1:s(1),ones(s(1),1),'r');
        v = axis;
        %axis([1,s(1),0,min(max(v(4),1.2),4)])
        axis tight
    end

    hi95=zeros(s(2),s(1)); % registro de faltas
    hi=zeros(s(2),s(1));
    for i=1:s(2),
        hi(i,:) = hotelling(:,i)>limb';
        hi95(i,:) = hotelling(:,i)>lima';
    end
else
    if opt,
        axes(axes1)
        hold off  
        plot(0,0)
    end
end
    
    

% SPE 

limar=[];
limbr=[];

% Estimation of the SPE control limits following Box's approximation
% approach
for i=1:s(1),   
   
   % Estimation of the SPE control limits following Box's approximation
    sumRes = sum(resmod(:,:,i).^2,2);    
    m = mean(sumRes);
    v = var(sumRes);
    %limar=[limar (v/(2*m))*chi2inv(1-alpr95,(2*m^2)/v)];         
    %limbr=[limbr (v/(2*m))*chi2inv(1-alpr,(2*m^2)/v)];             
   % Estimation of the SPE control limits following Jackson & Mudholkar's
   % approach
     limar=[limar spe_lim(resmod(:,:,i),0.05)]; % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
     limbr=[limbr spe_lim(resmod(:,:,i),0.01)]; % limite al 99% de conf.
end

res = residuals;

if opt,
    axes(axes2)
    hold off
    plot(res(1:s_sBn),tg)
    hold on;

    xlabel('Sampling time','FontSize', 12,'FontWeight','Bold');
    ylabel('SPE','FontSize', 12,'FontWeight','Bold');

    plot(1:s(1),real(limar),'r--','LineWidth',1.2);
    plot(1:s(1),real(limbr),'r','LineWidth',1.2);
    v = axis;
    v2=max(real(limbr));
    axis([1,s(1),0,min(max(v(4),v2),4*v2)])
end

ri95=zeros(s(2),s(1));
ri=zeros(s(2),s(1));
for i=1:s(2),
    ri(i,:) = residuals(:,i)>limbr';
    ri95(i,:) = residuals(:,i)>limar';
end

