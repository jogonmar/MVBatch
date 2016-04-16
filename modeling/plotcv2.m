function [alph,alpr,alph95,alpr95]=plotcv2(resmod,hotellingcv,residualscv,lotes,tg,pc,opt,axes1,axes2)

% Plots the D-statistic and SPE values of the calibration batches using 
%   leave-one-out cross-validation. 
%
% [alph,alpr,alph95,alpr95]=plotcv(resmod,hotellingcv,residualscv,lotes,tg
%   ,pc,opt) % call with standard parameters
%
% [alph,alpr,alph95,alpr95]=plotcv(resmod,hotellingcv,residualscv,lotes,tg
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
% tg: string with the plot options for the statistics of the test batch.
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
%
% codified by: José Camacho Páez and J.M. Gonzalez-Martinez.
% version: 0.0
% last modification: August/14.

% Parameters checking

if nargin < 6, error('Numero de argumentos erroneos.'); end;
if nargin < 7, opt = 1; end;
if nargin < 8, 
    h = figure;
    axes1 = axes; 
end;
if nargin < 9, 
    h2 = figure;
    axes2 = axes; 
end;

alph = 0.01;
alpr = 0.01;
alph95 = 0.05;
alpr95 = 0.05; 

% D-statistic

s=size(hotellingcv);
out = false;

lclu = s(1);

if pc(1)~=0,
    lima=[];
    limb=[];
    for i=1:s(1),
        lima(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-0.05,pc(i),lotes-pc(i)); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
        limb(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-0.01,pc(i),lotes-pc(i)); % limite al 99% de conf.
    end
    
    % experimental limit
    
    esc=[];
    for i=1:s(2),
        esc = [esc ; hotellingcv(:,i)./lima'];
    end
    esc2 = sort(esc);
                        
    alp = 0.05;
    lev = round(s(1)*s(2)*(1-alp));
    ind = esc2(lev);
    ind2 = mod(find(esc==ind,1),s(1));
    if ~ind2, ind2 = s(1); end;
    alph95=1-fcdf(((lotes*(lotes-pc(ind2)))/(pc(ind2)*(lotes*lotes-1)))*(lima(ind2) * ind),pc(ind2),lotes-pc(ind2));

    esc=[];
    for i=1:s(2),
        esc = [esc ; hotellingcv(:,i)./limb'];
    end
    esc2 = sort(esc);
    
    alp = 0.01;
    lev = round(s(1)*s(2)*(1-alp));
    ind = esc2(lev);
    ind2 = mod(find(esc==ind,1),s(1));
    if ~ind2, ind2 = s(1); end;
    
    alph=1-fcdf(((lotes*(lotes-pc(ind2)))/(pc(ind2)*(lotes*lotes-1)))*(limb(ind2) * ind),pc(ind2),lotes-pc(ind2));

    for i=1:s(1),
        limacv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alph95,pc(i),lotes-pc(i)); % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
        limbcv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alph,pc(i),lotes-pc(i)); % limite al 99% de conf.
    end

    if opt,
        axes(axes1)
        hold off
        plot(hotellingcv,tg);
        hold on

        xlabel('Sampling time','FontSize', 12,'FontWeight','bold');
        ylabel('D-statistic','FontSize', 12,'FontWeight','bold');

        plot(1:s(1),lima,'r--');
        plot(1:s(1),limb,'r');
        
        plot(1:s(1),limacv,'g--');
        plot(1:s(1),limbcv,'g');

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

limar=[];
limbr=[];

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

% experimental limit
    
    esc=[];
    for i=1:s(2),
        esc = [esc ; residualscv(:,i) ./ limar'];
    end
    esc2 = sort(esc);
                        
    alp = 0.05;
    lev = round(s(1)*s(2)*(1-alp));
    ind = esc2(lev);
    ind2 = mod(find(esc==ind,1),s(1));
    if ~ind2, ind2 = s(1); end;
    
    % Adjust of the confidence level to meet the imposed 95% confidence
    % level following Box's approximation
    %alpr95 = 1-chi2cdf((limar(ind2) * ind)/(v/(2*m)),(2*m^2)/v);
    
    % Adjust of the confidence level to meet the imposed 95% confidence
    % level following Jackson & Mudholkar's approach
    alpr95=spe_pvalue(resmod(:,:,ind2),(limar(ind2) * ind));
    
    esc=[];
    for i=1:s(2),
        esc = [esc ; residualscv(:,i) ./ limbr'];
    end
    esc2 = sort(esc);
                        
    alp = 0.01;
    lev = round(s(1)*s(2)*(1-alp));
    ind = esc2(lev);
    ind2 = mod(find(esc==ind,1),s(1));
    if ~ind2, ind2 = s(1); end;
    
    % Adjust of the confidence level to meet the imposed 95% confidence
    % level following Box's approximation
    %alpr = 1-chi2cdf((limbr(ind2) * ind)/(v/(2*m)),(2*m^2)/v);
    
     % Adjust of the confidence level to meet the imposed 99% confidence
    % level following Jackson & Mudholkar's approach
    alpr=spe_pvalue(resmod(:,:,ind2),(limbr(ind2) * ind));
  
    limarcv=[];
    limbrcv=[];

    % Estimation of the SPE control limits using CV follwing Box's
    % approximation
%     for i=1:s(1),
%         sumRes = sum(resmod(:,:,i).^2,2);    
%         m = mean(sumRes);
%         v = var(sumRes);
%         limarcv=[limarcv (v/(2*m))*chi2inv(1-alpr95,(2*m^2)/v)]; % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
%         limbrcv=[limbrcv (v/(2*m))*chi2inv(1-alpr,(2*m^2)/v)]; % limite al 99% de conf.
%     end 

      % Estimation of the SPE control limits using CV Jackson & Mudholkar's approach
     for i=1:s(1),
         limarcv=[limarcv spe_lim(resmod(:,:,i),alpr95)]; % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
         limbrcv=[limbrcv spe_lim(resmod(:,:,i),alpr)]; % limite al 99% de conf.
     end

if opt,
    axes(axes2)
    hold off
    plot(residualscv,tg);
    hold on;

    xlabel('Sampling time','FontSize', 12,'FontWeight','bold');
    ylabel('SPE','FontSize', 12,'FontWeight','bold');

    plot(real(limar),'r--');
    plot(real(limbr),'r');

    plot(real(limarcv),'g--');
    plot(real(limbrcv),'g');

    v=axis;
    axis([1 lclu v(3) v(4)]);
end