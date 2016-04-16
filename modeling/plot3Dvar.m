function [] = plot3Dvar(X,J,X2,varNames,savingpath)

% Funci???n que permite graficar las trayectorias de una determinada variable
% para todos los lotes.
if iscell(X), nbatches = length(X);
else
    nbatches = size(X,3);
end
   
plot2_1 = false; plot2_2 = false;
if nargin < 2, error('Number of argument should be two. '); end
if nargin < 3, X2 = []; end
if nargin > 2 && (iscell(X2) && ~isempty(X2)), plot2_1 = true;
elseif nargin > 2 && (~iscell(X2) && ~isempty(X2)),plot2_2 = true;
end
if nargin <4, varNames=cell(0,0);end
saving = 1;
if nargin <5, saving = 0;end

fig_h = figure;

if iscell(X)

   %figure;
   % grey [0.745098 0.745098 0.745098]
   for i=1:nbatches
      plot(X{i}(:,J),'-','Color',[0.941176 0.972549 1],'LineWidth',1); hold on;
   end
   if plot2_1
       for i=1:length(X2)
            plot(X2{i}(:,J),'k','LineWidth',1);
       end
   elseif plot2_2
        plot(X2(:,J),'k','LineWidth',1.5);
      
   end
else
    %figure;
    var = squeeze(X(:,J,:));

    for i = 1:size(var,2)
        plot(var(:,i),'-','Color',[0.941176 0.972549 1],'LineWidth',1); hold on;
    end
       if plot2_2
           var2 = squeeze(X2(:,J,:));
           for i=1:size(X2,3)
                plot(var2(:,i),'k-','LineWidth',1.5);
           end
       end
end

xlabel('Batch time','FontSize',14);
if ~isempty(varNames)
    ylabel(varNames{J,1},'FontSize',14);
else
    ylabel('Process variable','FontSize',14);
end
axes_h=get(fig_h,'Children');
set(axes_h(1),'FontSize',12);
set(fig_h,'Color','w');
axis tight    

if saving
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2',savingpath,'-loose')
end


