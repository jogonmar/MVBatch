function  [varargout] = plot3D_batchtools(x,test,vars,varNames,uipanelPlots)

% Plots 3-way batch data for each process variable in the layout of the GUI. 
%
% CALLS:
%
% hnd = plot3D_batchtools(x,test,vars)                        % minimum call
% hnd = plot3D_batchtools(x,test,vars,varNames,uipanelPlots)  % complete call 
%
% INPUTS:
%
% x: (Ix1) cell structure containing the field data, which contains z cell array with batch data sampled at certain frequency
%
% test: batch to highlight in the visualization with the same structure as
% x.
%
% vars: (nx1) vector with the indices of the process variables to plot. 
%
% varNames: cell array with the tagnames of the process variables.
%
% uipanelPlots: (1x1) handle to the figure.
%
%
% OUTPUTS:
%
% varargout: (1x1) handle to the figure.
%
% coded by: José M. González-Martínez (jogonmar@gmail.com)       
% last modification: Aug/14.
%
% Copyright (C) 2016  Technical University of Valencia, Valencia
% Copyright (C) 2016  José M. González-Martínez
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
if nargin < 3, error('Error in the number of arguments.'); end;
if ~isfield(x{1,1},'data'), error('Incorrect field in the data structure.'); end
if ~isfield(test{1,1},'data'), error('Incorrect field in the data structure.'); end
if ~isvector(vars), error('Incorrect type of data structure for vars.');end
sv = numel(vars);
if nargin < 4, varNames = repmat({'Variable ','Sampling point','Units'},sv,1); end
if nargin < 5, uipanelPlots = figure; end

if sv >4 
    columns = 3;
    rows = 3;
elseif sv >=3 && sv <=4
    columns = 2;
    rows = 2;
elseif sv ==2
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
    for i=1:length(x),
        for z=1:size(x{i}.data,1)
            maxy = nanmax(x{1}.data{z}(:,indj(j)+2));
            miny = nanmin(x{1}.data{z}(:,indj(j)+2));
            maxk = Inf;
            plot(x{i}.data{z}(:,1),x{i}.data{z}(:,indj(j)+2),'b-','Color',[0.466667 0.533333 0.68]); hold on;
            plot(test{1}.data{z}(:,1),test{1}.data{z}(:,indj(j)+2),'r-','LineWidth',1);
            maxy = nanmax(nanmax(x{1}.data{z}(:,indj(j)+2)),maxy);
            miny = nanmin(nanmin(x{1}.data{z}(:,indj(j)+2)),miny);
            maxk = nanmax(maxk,size(x{1}.data{z},1));
        end
    end
    title(varNames(vars(indj(j)),1),'fontweight','b');
    xlabel(varNames(vars(indj(j)),2));
    ylabel(varNames(vars(indj(j)),3));
    if miny==maxy, axis([1 maxk miny-1 maxy+1]);
    else
        if ~isnan(miny) && ~isnan(maxy), axis([1 maxk miny maxy]); end
    end
axis tight
end

if nargout>0,   
    varargout{1} = hnd;
end
