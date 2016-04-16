function [alph,alpr,alph95,alpr95]=plotOverall(resmod,hotelling,residuals,lotes,pc,opt,axes1,axes2)

% Plots the overall D-statistic and SPE values of the calibration batches using 
%   the leave-one-out cross-validated and theorical control limits. 
%
% [alph,alpr,alph95,alpr95]=plotcv(resmod,hotelling,residuals,lotes,tg
%   ,pc,opt) % call with standard parameters
%
% [alph,alpr,alph95,alpr95]=plotcv(resmod,hotelling,residuals,lotes,tg
%   ,pc,opt,axes1,axes2) % complete call
%
%
% INPUTS:
%
% resmod: (IxJxK) residuals in the calibration data, K(sampling times) 
%       x J(variables) x I(batches)
%
% hotelling_cv: (Kx1) D-statistic of the calibration batches obtained in a
%    leave-one-out cross-validation, K(sampling times).
%
% residuals_cv: (Kx1) SPE of the calibration batches obtained in a
%    leave-one-out cross-validation, K(sampling times).
%
% lotes: (1x1) number of calibration batches.
%
% pc: (1x1) number of principal components for the D-statistic.
%
% opt: boolean (1x1) 
%       true: plot results.
%       false: do not plot results.
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
% codified by: José M. Gonzalez-Martinez.
% version: 0.0
% last modification: 20/Oct/13.

% Parameters checking

if nargin < 5, error('The number of argument is not correct.'); end;
if nargin < 6, opt = 1; end;
if nargin < 7, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 8, 
    h2 = figure;
    axes2 = axes; 
end;

alph = 0.01;
alpr = 0.01;
alph95 = 0.05;
alpr95 = 0.05; 

% D-statistic

s=size(hotelling);

lclu = s(1);

if pc(1)~=0,
    lima =(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-0.05,pc,lotes-pc); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
    limb =(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-0.01,pc,lotes-pc); % limite al 99% de conf.
    
    % experimental limit
    
    esc=[];
    esc = [esc ; hotelling./lima'];
    esc2 = sort(esc);
                        
    alp = 0.05;
    lev = round(s(1)*(1-alp));
    ind = esc2(lev);
    ind2 = lev;
    alph95=1-fcdf(((lotes*(lotes-pc))/(pc*(lotes*lotes-1)))*(lima * ind),pc,lotes-pc);

    esc=[];
    esc = [esc ; hotelling./limb'];
    esc2 = sort(esc);
    
    alp = 0.01;
    lev = round(s(1)*(1-alp));
    ind = esc2(lev);
    ind2 = lev;
    alph=1-fcdf(((lotes*(lotes-pc))/(pc*(lotes*lotes-1)))*(limb * ind),pc,lotes-pc);

 
    limacv=(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-alph95,pc,lotes-pc); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
    limbcv=(pc*(lotes*lotes-1)/(lotes*(lotes-pc)))*finv(1-alph,pc,lotes-pc); % limite al 99% de conf. 
   
    if opt,
        axes(axes1)
        hold off
        bar(hotelling,'b');
        hold on

        xlabel('Batches','FontSize', 12,'FontWeight','bold');
        ylabel('D-statistic','FontSize', 12,'FontWeight','bold');

        plot(1:s(1),repmat(lima,s(1),1),'r--');
        plot(1:s(1),repmat(limb,s(1),1),'r');
        
        plot(1:s(1),repmat(limacv,s(1),1),'g--');
        plot(1:s(1),repmat(limbcv,s(1),1),'g');

        v=axis;
        axis([1 lclu v(3) v(4)]);
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
    %limar= spe_lim(resmod,0.05); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
    %limbr= spe_lim(resmod,0.01); % limite al 99% de conf.  
    
   % Estimation of the SPE control limits following Box's approximation
   % approach
    E = unfold(permute(resmod,[3,2,1]),Inf);
    m = mean(sum(E.^2,2));
    v = var(sum(E.^2,2));
    limar = (v/(2*m))*chi2inv(1-alpr95,(2*m^2)/v);
    limbr = (v/(2*m))*chi2inv(1-alpr,(2*m^2)/v);
    
% experimental limit
    
    esc=[];
    for i=1:s(2),
        esc = [esc ; residuals./ limar'];
    end
    esc2 = sort(esc);
                        
    alp = 0.05;
    lev = round(s(1)*(1-alp));
    ind = esc2(lev);
    ind2 = lev;
    
    % Adjust of the confidence level to meet the imposed 95% confidence
    % level following Box's approximation
    alpr95 = 1-chi2cdf((limar*ind)/(v/(2*m)),(2*m^2)/v);
    
    % Adjust of the confidence level to meet the imposed 95% confidence
    % level following Jackson & Mudholkar's approach
    %alpr95=spe_pvalue(resmod,(limar * ind));
    
    esc=[];
    for i=1:s(2),
        esc = [esc ; residuals ./ limbr'];
    end
    esc2 = sort(esc);
                        
    alp = 0.01;
    lev = round(s(1)*(1-alp));
    ind = esc2(lev);
    ind2 = lev;    
    
    % Adjust of the confidence level to meet the imposed 99% confidence
    % level following Box's approximation
    alpr = 1-chi2cdf((limbr*ind)/(v/(2*m)),(2*m^2)/v);
    
    % Adjust of the confidence level to meet the imposed 99% confidence
    % level following Jackson & Mudholkar's approach
    
%    alpr=spe_pvalue(resmod,(limbr * ind));
   
    % Estimation of the SPE control limits using CV follwing Box's
    % approximation
    
     limarcv = (v/(2*m))*chi2inv(1-alpr95,(2*m^2)/v);
     limbrcv = (v/(2*m))*chi2inv(1-alpr,(2*m^2)/v);
    
    % Estimation of the SPE control limits using CV Jackson & Mudholkar's approach
    %limarcv= spe_lim(unfold(permute(resmod,[3 2 1]),Inf),alpr95); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
    %limbrcv= spe_lim(unfold(permute(resmod,[3 2 1]),Inf),alpr); % limite al 99% de conf.  

if opt,
    axes(axes2)
    hold off
    bar(residuals,'b');
    hold on;
    xlabel('Batches','FontSize', 12,'FontWeight','bold');
    ylabel('SPE','FontSize', 12,'FontWeight','bold');
    
    plot(1:s(1),repmat(real(limar),s(1),1),'r--');
    plot(1:s(1),repmat(real(limbr),s(1),1),'r');

    plot(1:s(1),repmat(real(limarcv),s(1),1),'g--');
    plot(1:s(1),repmat(real(limbrcv),s(1),1),'g');

    v=axis;
    axis([1 lclu v(3) v(4)]);
end