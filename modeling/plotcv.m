function [alph,alpr,alph95,alpr95]=plotcv(resmod,hotellingcv,residualscv,lotes,tg,pc,opt,axes1,axes2)

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
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 19/May/09.

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
    hi = length(find(hotellingcv>lima'*ones(1,s(2))));
    alp = 0.05;
    lev = round(s(1)*s(2)*alp);
    sal = false;
    c=0;
    paso=0.001;
    for i=1:s(1),
        limacv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alp,pc(i),lotes-pc(i));
    end
    while hi~=lev && ~sal,
        if hi > lev,
            while alp<=paso,
                paso = 0.1*paso;
            end
            alp = alp - paso;
            if c==-1,
                sal=true;
            end
            c=1;
        else
            while alp>=1-paso,
                paso = 0.1*paso;
            end
            alp = alp + paso;
            if c==1,
                sal=true;
            end
            c=-1;
        end
        for i=1:s(1),
            limacv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alp,pc(i),lotes-pc(i));
        end
        hi = length(find(hotellingcv>limacv'*ones(1,s(2))));
    end
    alph95=alp;

    hi = length(find(hotellingcv>limb'*ones(1,s(2))));
    alp = 0.01;
    lev = round(s(1)*s(2)*alp);
    sal = false;
    c=0;
    paso=0.001;
    for i=1:s(1),
        limbcv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alp,pc(i),lotes-pc(i));
    end
    while hi~=lev && ~sal,
        if hi > lev,
            while alp<=paso,
                paso = 0.1*paso;
            end
            alp = alp - paso;   
            if c==-1,
                sal=true;
            end
            c=1;
        else
            while alp>=1-paso,
                paso = 0.1*paso;
            end
            alp = alp + paso;
            if c==1,
                sal=true;
            end
            c=-1;
        end
        for i=1:s(1),
            limbcv(i)=(pc(i)*(lotes*lotes-1)/(lotes*(lotes-pc(i))))*finv(1-alp,pc(i),lotes-pc(i));
        end
        hi = length(find(hotellingcv>limbcv'*ones(1,s(2))));
    end
    
    if opt,
        axes(axes1)
        hold off
        plot(hotellingcv,tg);
        hold on

        xlabel('Sampling time','FontSize', 14);
        ylabel('D-statistic','FontSize', 14);

        plot(1:s(1),lima,'r--');
        plot(1:s(1),limb,'r');
        
        plot(1:s(1),limacv,'g--');
        plot(1:s(1),limbcv,'g');

        v=axis;
        axis([1 lclu v(3) v(4)]);
    end
    alph=alp;
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
    limar=[limar spe_lim(resmod(:,:,i),0.05)]; % limite al 95% de conf., importante pasarle a residuallimit una matriz para que identifique residuos
    limbr=[limbr spe_lim(resmod(:,:,i),0.01)]; % limite al 99% de conf.
end

% experimental limit

limarcv=limar;
hi = length(find(residualscv>limar'*ones(1,s(2))));
alp = 0.05;
lev = round(s(1)*s(2)*alp);
sal = false;
c=0;
paso=0.001;
while hi~=lev && ~sal,
    if hi > lev,
        while alp<=paso,
            paso = 0.1*paso;
        end
        alp = alp - paso; 
        if c==-1,
            sal=true;
        end
        c=1;
    else
        while alp>=1-paso,
            paso = 0.1*paso;
        end
        alp = alp + paso;
        if c==1,
            sal=true;
        end
        c=-1;
    end

    for i=1:s(1),
        limarcv(i)=spe_lim(resmod(:,:,i),alp); 
    end
    hi = length(find(residualscv>limarcv'*ones(1,s(2))));
end
alpr95=alp;

limbrcv=limar;
hi = length(find(residualscv>limbr'*ones(1,s(2))));
alp = 0.01;
lev = round(s(1)*s(2)*alp);
sal = false;
c=0;
paso=0.001;
while hi~=lev && ~sal,
    if hi > lev,
        while alp<=paso,
            paso = 0.1*paso;
        end
        alp = alp - paso;   
        if c==-1,
            sal=true;
        end
        c=1;
    else
        while alp>=1-paso,
            paso = 0.1*paso;
        end
        alp = alp + paso;
        if c==1,
            sal=true;
        end
        c=-1;
    end

    for i=1:s(1),
        limbrcv(i)=spe_lim(resmod(:,:,i),alp); 
    end
    hi = length(find(residualscv>limbrcv'*ones(1,s(2))));
end
alpr=alp;

if opt,
    axes(axes2)
    hold off
    plot(residualscv,tg);
    hold on;

    xlabel('Sampling time','FontSize', 14);
    ylabel('SPE','FontSize', 14);

    plot(real(limar),'r--');
    plot(real(limbr),'r');

    plot(real(limarcv),'g--');
    plot(real(limbrcv),'g');

    v=axis;
    axis([1 lclu v(3) v(4)]);
end