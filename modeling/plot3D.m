function  [varargout] = plot3D(x,clu,test,rows,vars)

% Plots 3-way batch data. 
%
% plot3D(x)           % plot batch data
% plot3D(x,clu,test)       % plot batch data, phases and a test batches
% hnd = plot3D(x,clu,test,rows,vars)  % complete call 
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
% test: can be one of two possibilities (consistent with x):
%
%   - (KxJxIt) three-way batch data matrix, K(sampling times) x J(variables)
%       x It(batches)
%
%   - cell (1xIt) group of batches with a common set of variables (some 
%       of the variables may not have been measured in some batches). In 
%       this case, x{i} (KixJ2) contains the data corresponding to the J2 
%       variables collected in the batch i with a common sampling timing. 
%
%
% rows: number of rows in the image (3 by default);
%
% vars: index to the vars to plot (all by default);
%
%
% OUTPUTS:
%
% hnd: (1x1) handle to the figure/s.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/Aug/09
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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

if nargin < 1, error('Error in the number of arguments.'); end;
s = size(x);

if ~iscell(x),
    if ~isstruct(x),
        data_type = 0;
        if length(s)<3, s(3)=1; end;
    else
        data_type = 1;
    end
else
    data_type = 2;
end
    
if nargin < 2 || isempty(clu), 
    clu=1;
end

if nargin < 3 || isempty(test), 
    tt = false; 
else
    tt = true;
end;

if nargin < 4, rows = 3; end;
columns = round(rows*4/3);
max_var = rows*columns;

if nargin == 5, 
    switch data_type,
        case 1,
            s2 = 0;
            ind_v=[];
            ind_v2=[];
            for i=1:length(x(1).data),
                s2=s2+size(x(1).data{i},2)-1;
                ind_v=[ind_v i*ones(1,size(x(1).data{i},2)-1)];
                ind_v2=[ind_v2 2:(size(x(1).data{i},2))];
            end
            for i=1:length(x),
                j_new=1;
                j=vars(1);
                while j<=vars(end),
                    if find(j==vars),
                        ind_j = find(ind_v(j)==ind_v(vars));
                        y(i).data{j_new}=x(i).data{ind_v(j)}(:,[1 ind_v2(vars(ind_j))]);
                        j_new = j_new+1;
                        j = vars(ind_j(end));
                    end
                    j = j+1;
                end
            end
            x=y;
            
            if tt,            
                for i=1:length(test),
                    j_new=1;
                    j=vars(1);
                    while j<=vars(end),
                        if find(j==vars),
                            ind_j = find(ind_v(j)==ind_v(vars));
                            ty(i).data{j_new}=test(i).data{ind_v(j)}(:,[1 ind_v2(vars(ind_j))]);
                            j_new = j_new+1;
                            j = vars(ind_j(end));
                        end
                        j = j+1;
                    end
                end
                test=ty;
            end
                

        case 2,
            for i=1:length(x),
                x{i} = x{i}(:,vars);
            end
            
            if tt,
                for i=1:length(test),
                    test{i} = test{i}(:,vars);
                end
            end

        otherwise
            x = x(:,vars,:);
            
            if tt,
                test = test(:,vars,:);
            end
    end;
    
    s = size(x);
    if data_type == 0 && length(s)<3, s(3)=1; end;
    
end;

% Computation

hnd=[];

switch data_type,
    case 1,
        s2 = 0;
        ind_v=[];
        ind_v2=[];
        s1 = [];
        for i=1:length(x(1).data),
            s2=s2+size(x(1).data{i},2)-1;
            ind_v=[ind_v i*ones(1,size(x(1).data{i},2)-1)];
            ind_v2=[ind_v2 2:(size(x(1).data{i},2))];
            s1=[s1 x(1).data{i}(end,1)*ones(1,size(x(1).data{i},2)-1)];
        end
        maxk = s1;
        if max_var<maxk,
            columns = ceil(max_var/rows);
        end

        for k=1:max_var:s2,
            hnd = [hnd figure];
            indj = k:nanmin(s2,k+max_var-1);
            lenj = length(indj);
            for j=1:lenj,
                subplot(rows,columns,j);
                maxy = nanmax(x(1).data{ind_v(indj(j))}(:,ind_v2(indj(j))));
                miny = nanmin(x(1).data{ind_v(indj(j))}(:,ind_v2(indj(j)))); 
                hold on
                for i=1:nanmax(s),
                    plot(x(i).data{ind_v(indj(j))}(:,1),x(i).data{ind_v(indj(j))}(:,ind_v2(indj(j))),'b-')
                    maxk(j) = nanmax(maxk(j),x(i).data{ind_v(indj(j))}(end,1));
                    maxy = nanmax(maxy,nanmax(x(i).data{ind_v(indj(j))}(:,ind_v2(indj(j)))));
                    miny = nanmin(miny,nanmin(x(i).data{ind_v(indj(j))}(:,ind_v2(indj(j)))));
                end
                if tt,
                    st = size(test);                   
                    for i=1:nanmax(st),
                        plot(test(i).data{ind_v(indj(j))}(:,1),test(i).data{ind_v(indj(j))}(:,ind_v2(indj(j))),'r-','LineWidth',1.5)
                        maxk(j) = nanmax(maxk(j),test(i).data{ind_v(indj(j))}(end,1));
                        maxy = nanmax(maxy,nanmax(test(i).data{ind_v(indj(j))}(:,ind_v2(indj(j)))));
                        miny = nanmin(miny,nanmin(test(i).data{ind_v(indj(j))}(:,ind_v2(indj(j)))));
                    end
                end
                if maxy==miny, maxy=maxy+0.5; miny=miny-0.5; end;
                for i=2:nanmax(clu),
                    ind=find(clu==i,1);
                    plot([ind ind],[maxy miny],'k--');
                end
                if maxk(j)>1, axis([1 maxk(j) miny maxy]); end;
                axes_h=get(hnd(end),'Children');
                ticks=get(axes_h,'XTick');
                tickl=get(axes_h,'XTickLabel');
                if iscell(ticks),
                    set(axes_h,'XTick',ticks{1}([1:rows:end-rows end]))
                    set(axes_h,'XTickLabel',tickl{1}([1:rows:end-rows end],:))
                else
                    set(axes_h,'XTick',ticks([1:rows:end-rows end]))
                    set(axes_h,'XTickLabel',tickl([1:rows:end-rows end],:))
                end
                axis tight
            end
        end
     
    case 2,
        s2=size(x{1});
        maxk = s2(1);
        if max_var>s2(2),
            columns = ceil(max_var/rows);
        end

        for k=1:max_var:s2(2),
            hnd = [hnd figure];
            indj = k:nanmin(s2(2),k+max_var-1);
            lenj = length(indj);
            for j=1:lenj,
                subplot(rows,columns,j);
                maxy = nanmax(x{1}(:,indj(j)));
                miny = nanmin(x{1}(:,indj(j)));
                hold on
                for i=1:nanmax(s),
                    plot(x{i}(:,indj(j)),'b-')
                    maxk = nanmax(maxk,size(x{i},1));
                    maxy = nanmax(maxy,nanmax(x{i}(:,indj(j))));
                    miny = nanmin(miny,nanmin(x{i}(:,indj(j))));
                end
                if tt,
                    st = size(test);
                    for i=1:nanmax(st),
                        plot(test{i}(:,indj(j)),'r-','LineWidth',1.5)
                        maxk = nanmax(maxk,size(test{i},1));
                        maxy = nanmax(maxy,nanmax(test{i}(:,indj(j))));
                        miny = nanmin(miny,nanmin(test{i}(:,indj(j))));
                    end
                end
                if maxy==miny, maxy=maxy+0.5; miny=miny-0.5; end;
                for i=2:nanmax(clu),
                    ind=find(clu==i,1);
                    plot([ind ind],[maxy miny],'k--');
                end
               axis([1 maxk miny maxy]);
                axes_h=get(hnd(end),'Children');
%                 ticks=get(axes_h,'XTick');
%                 tickl=get(axes_h,'XTickLabel');
%                 if iscell(ticks),
%                     set(axes_h,'XTick',ticks{1}([1:rows:end-rows end]))
%                     set(axes_h,'XTickLabel',tickl{1}([1:rows:end-rows end],:))
%                 else
%                     set(axes_h,'XTick',ticks([1:rows:end-rows end]))
%                     set(axes_h,'XTickLabel',tickl([1:rows:end-rows end],:))
%                 end
                axis tight
            end
        end
        
    otherwise,
        if max_var>s(2),
            columns = ceil(max_var/rows);
        end
        
        for k=1:max_var:s(2),
            hnd = [hnd figure];
            indj = k:nanmin(s(2),k+max_var-1);
            lenj = length(indj);
            for j=1:lenj,
                subplot(rows,columns,j);
                m = nanmin(nanmin(x(:,indj(j),:)));
                M = nanmax(nanmax(x(:,indj(j),:)));
                hold on
                for i=1:s(3),
                    plot(x(:,indj(j),i),'b');
                end
                if tt,
                    st = size(test);
                    if ndims(st)<3, st=[st 1]; end;
                    for i=1:st(3),
                        plot(test(:,indj(j),i),'r-')
                    end
                    m = nanmin(m,nanmin(nanmin(test(:,indj(j),:))));
                    M = nanmax(M,nanmax(nanmax(test(:,indj(j),:))));
                end
                
                if M==m, M=M+0.5; m=m-0.5; end;
                for i=2:nanmax(clu),
                    ind=find(clu==i);
                    plot([ind(1) ind(1)],[M m],'k--');
                end
%                axis([1 s(1) m M]);
%                 axes_h=get(hnd(end),'Children');
%                 ticks=get(axes_h,'XTick');
%                 tickl=get(axes_h,'XTickLabel');
%                 if iscell(ticks),
%                     set(axes_h,'XTick',ticks{1}([1:rows:end-rows end]))
%                     set(axes_h,'XTickLabel',tickl{1}([1:rows:end-rows end],:))
%                 else
%                     set(axes_h,'XTick',ticks([1:rows:end-rows end]))
%                     set(axes_h,'XTickLabel',tickl([1:rows:end-rows end],:))
%                 end
                axis tight
            end
        end
end

if nargout>0,   
    varargout{1} = hnd;
end

