function  [varargout] = plot3D_batchtools(x,clu,obs,vars,varNames,test, uipanelPlots)

% Plots 3-way batch data. 
%
% plot3D(x)           % plot batch data
% plot3D(x,clu,handles.test)       % plot batch data, phases and a handles.test batches
% hnd = plot3D(x,clu,handles.test,rows,vars)  % complete call 
%
% INPUTS:
%
% x: can be one of two possibilities:
%
%   - (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%       x I(batches)
%
%   - cell (1xI) group of batches with a common set of variables (some 
%       of the variables may not have been measured in some batches). In 
%       this case, x{i} (KixJ2) contains the data corresponding to the J2 
%       variables collected in the batch i with a common sampling timing. 
%
% clu: (1xK) vector with the assignment of the sampling times to the phases, 
%   numbered from 1 onwards (1 phase by default).
%
% obs: (nx1) vector with the indices of the batches to plot. 
%
% vars: (nx1) vector with the indices of the process variables to plot. 
%
% varNames: cell array with the tagnames of the process variables.
%
% test: batch to highlight in the visualization
%
% uipanelPlots: (1x1) handle to the figure.
%
%
% OUTPUTS:
%
% varargout: (1x1) handle to the figure.
%
% version: 1.0
% last modification: Aug/14.

% Parameters checking
if nargin < 1, error('Error in the number of arguments.'); end;
s = size(obs);
    
if nargin < 2 || isempty(clu), 
    clu=1;
end
if nargin < 3, vars = 9;end

%if nargin < 4, rows = 3; end;
if numel(vars) >4 
    columns = 3;
    rows = 3;
elseif numel(vars) >=3 && numel(vars) <=4
    columns = 2;
    rows = 2;
elseif numel(vars) ==2
    columns = 1;
    rows = 2;
else
    columns = 1;
    rows = 1;
end
    
 max_var = rows*columns;

% Computation


s2=size(vars);

indj = 1:min(s2(1),max_var);
lenj = length(indj);
hnd = zeros(lenj,1);

for j=1:lenj,
    hnd(j) = subplot(rows,columns,j, 'Parent',uipanelPlots);
    cla reset
    hold on
    for i=1:max(s),
        for z=1:size(x{i}.data,1)
            maxy = nanmax(x{1}.data{z}(:,indj(j)+1));
            miny = nanmin(x{1}.data{z}(:,indj(j)+1));
            maxk = Inf;
            plot(x{i}.data{z}(:,1),x{i}.data{z}(:,indj(j)+1),'b-','Color',[0.466667 0.533333 0.68]); hold on;
            plot(test{1}.data{z}(:,1),test{1}.data{z}(:,indj(j)+1),'r-','LineWidth',1);

            maxy = nanmax(nanmax(x{1}.data{z}(:,indj(j)+1)),maxy);
            miny = nanmin(nanmin(x{1}.data{z}(:,indj(j)+1)),miny);
            maxk = nanmax(maxk,size(x{1}.data{z},1));
            
        end
    end
    title(varNames(vars(indj(j)),1),'fontweight','b');
    xlabel(varNames(vars(indj(j)),2));
    ylabel(varNames(vars(indj(j)),3));
    for i=2:max(clu),
        ind=find(clu==i,1);
        plot([ind ind],[maxy miny],'k--');
    end
    if miny==maxy, axis([1 maxk miny-1 maxy+1]);
    else
        if ~isnan(miny) && ~isnan(maxy), axis([1 maxk miny maxy]); end
    end
axis tight
end

if nargout>0,   
    varargout{1} = hnd;
end
