function [nmap,map] = cov_map(x,x_ref,c_max,type,p_min,opt,console)

% Computes and plots batch-wise normalized covariance map.
%
% map = cov_map(x)          % total covariance map
% map = cov_map(x,x_ref,c_max,type,p_min,opt) % output in MATLAB console 
% map = cov_map(x,x_ref,c_max,type,p_min,opt,console)  % complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% x_ref: (KxJxI) reference data for normalization, K(sampling times) x 
%   J(variables) x I(batches)
%
% c_max: (1x1) value in the interval (0,1] to set the maximum value 
%   in the colormap (0.9 by default).
%
% type: (1x1) type of covariance map.
%       1: total covariance map (by default).
%       2: dynamic partial covariance map (pseudoinverse). 
%       3: instantaneous+dynamic partial covariance map (pseudoinverse). 
%
% p_min: (1x1) value in the interval (0,1] to set the minimum portion of
%   the total variance for an eigenvalue to be taken into account in types 
%   2 and 3 (0.1 by default).
%
% opt: boolean (1x1) 
%       true: plot results.
%       false: do not plot results.
%
% console: (1x1) handle of the EditText of the interface, 0 stands for the
%   MATLAB console (by default)
%
%
% OUTPUTS:
%
% nmap: (J·KxJ·K) normalized covariance map.
%
% map: (J·KxJ·K) covariance map.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 20/May/09
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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
if nargin < 2 || isempty(x_ref), x_ref = x; end; 
if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 3, c_max = 0.9; end;
if (c_max<=0||c_max>1), error('Incorrect value of c_max.'); end;
if nargin < 4, type = 1; end;
if (type<1||type>4), error('Incorrect value of type.'); end;
if nargin < 5, p_min = 0.1; end;
if (p_min<0||p_min>1), error('Incorrect value of p_min.'); end;
if nargin < 6, opt = 1; end;
if nargin < 7, console = 0; end;

% Computation

xu = preprocess2D(unfold(x,Inf),1);
xu_ref = preprocess2D(unfold(x_ref,Inf),1);
s = size(xu);
so = size(x);
       
cprintMV(console,'Processing.... Please, be patient.',[],0);
               
switch type,
    
    case 1, % Total Covariance Map
        map = (xu'*xu)/(s(1)-1);
        tit = 'Covariance Map';
 
    case 2, % MV-wise Partial (Backwards) Covariance Map using the Pseudo-inverse adjusting both for dynamic relationships
        tit = 'Dynamic Partial (Backwards) Covariance Map';
        map = (xu'*xu)/(s(1)-1);
        eu = xu; % backward regression
        au = xu; % fordward regression
            
        for i=1:so(1)-1, % time distance in MV units
            % Backward regression
            for o=s(2):-so(2):(i+1)*so(2),
                indk = o-(i+1)*so(2)+1:o-i*so(2);
                indk1 = o-so(2)+1:o;
                p = au(:,indk)'*au(:,indk); % Pseudo-inverse with specific tolerance
                est_eu(:,indk1) = au(:,indk)*(pinv(p,p_min*trace(p))*(au(:,indk)'*eu(:,indk1)));
                Ckk_eu(1:so(2),indk) = au(:,indk)'*eu(:,indk1)/(s(1)-1);  
            end
            
            % Map assignment
            diag_i = [];
            diag_s = [];    
            for o=1:so(1)-i, % measurement vector
                for u=1:so(2), % columns
                    ini_i = (so(2)*i+1) + ((o-1)*(s(2)+1)*so(2)) + (u-1)*s(2); 
                    ini_s = (s(2)*so(2)*i+1) + ((o-1)*(s(2)+1)*so(2)) + (u-1)*s(2); 
                    diag_i = [diag_i ini_i:ini_i+so(2)-1];
                    diag_s = [diag_s ini_s:ini_s+so(2)-1];
                end
            end
            Ckk_eu = Ckk_eu(1:so(2),1:s(2)-i*so(2));
            Ckk_eu = Ckk_eu(1:end);
            map(diag_s) = Ckk_eu;
            map(diag_i) = Ckk_eu;

            % Forward regression
            for o=1:so(2):s(2)-(i+1)*so(2),
                indk = o+i*so(2):o+so(2)*(1+i)-1;
                indk1 = o:o+so(2)-1;
                p = eu(:,indk)'*eu(:,indk); % Pseudo-inverse with specific tolerance
                au(:,indk1) = au(:,indk1) - eu(:,indk)*(pinv(p,p_min*trace(p))*(eu(:,indk)'*au(:,indk1)));
            end
            
            eu(:,(i+1)*so(2):s(2)) = eu(:,(i+1)*so(2):s(2)) - est_eu(:,(i+1)*so(2):s(2));
        end

    case 3, % MV-wise Partial (Backwards) Covariance Map using the Pseudo-inverse adjusting both for dynamic and instantaneous relationships
        tit = 'Partial (Backwards) Covariance Map';
        map = (xu'*xu)/(s(1)-1);
        eu = xu; % backward regression
        au = xu; % fordward regression
        
        for o=1:so(2):s(2),
            for i=0:so(2)-1, % instantaneous regression
                ind_rest = [o:o+i-1 o+i+1:o+so(2)-1];
                ind = o+i;
                p = au(:,ind_rest)'*au(:,ind_rest); % Pseudo-inverse with specific tolerance
                ini_eu(:,ind) = eu(:,ind) - au(:,ind_rest)*(pinv(p,p_min*trace(p))*(au(:,ind_rest)'*eu(:,ind)));
            end  
        end
            
        for i=1:so(1)-1, % time distance 

            % Forward regression without specific variable
            for o=1:so(2):s(2)-i*so(2),
                indk = o+i*so(2):o+so(2)*(1+i)-1;
                indk1 = o:o+so(2)-1;
                
                for j=0:so(2)-1, % instantaneous regression
                    ind_rest = indk([1:j j+2:end]);
                    p = eu(:,ind_rest)'*eu(:,ind_rest); % Pseudo-inverse with specific tolerance
                    est_au(:,indk1) = eu(:,ind_rest)*(pinv(p,p_min*trace(p))*(eu(:,ind_rest)'*au(:,indk1)));
                    ausp(:,indk1,j+1) = au(:,indk1) - est_au(:,indk1);
                end
                
                p = eu(:,indk)'*eu(:,indk); % Pseudo-inverse with specific tolerance
                est_au(:,indk1) = eu(:,indk)*(pinv(p,p_min*trace(p))*(eu(:,indk)'*au(:,indk1)));
            end          

            % Backward regression without specific variable
            for o=s(2):-so(2):(i+1)*so(2),
                indk = o-(i+1)*so(2)+1:o-i*so(2);
                indk1 = o-so(2)+1:o;
                
                for j=0:so(2)-1, % instantaneous regression
                    Ckk_eu(1:so(2),indk(1+j)) = ausp(:,indk,j+1)'*ini_eu(:,indk1(1+j))/(s(1)-1);  
                end
  
                p = au(:,indk)'*au(:,indk); % Pseudo-inverse with specific tolerance
                est_eu(:,indk1) = au(:,indk)*(pinv(p,p_min*trace(p))*(au(:,indk)'*eu(:,indk1)));
                eu(:,indk1) = eu(:,indk1) - est_eu(:,indk1); 
                for j=0:so(2)-1, % instantaneous regression
                    p = ausp(:,indk,j+1)'*ausp(:,indk,j+1); % Pseudo-inverse with specific tolerance
                    est_eu(:,indk1(1+j)) = ausp(:,indk,j+1)*(pinv(p,p_min*trace(p))*(ausp(:,indk,j+1)'*ini_eu(:,indk1(1+j))));
                end
                ini_eu(:,indk1) = ini_eu(:,indk1) - est_eu(:,indk1);
            end
            
            au(:,1:s(2)-i*so(2)) = au(:,1:s(2)-i*so(2)) - est_au(:,1:s(2)-i*so(2));
            
            % Map assignment
            diag_i = [];
            diag_s = [];    
            for o=1:so(1)-i, % measurement vector
                for u=1:so(2), % columns
                    ini_i = (so(2)*i+1) + ((o-1)*(s(2)+1)*so(2)) + (u-1)*s(2); 
                    ini_s = (s(2)*so(2)*i+1) + ((o-1)*(s(2)+1)*so(2)) + (u-1)*s(2); 
                    diag_i = [diag_i ini_i:ini_i+so(2)-1];
                    diag_s = [diag_s ini_s:ini_s+so(2)-1];
                end
            end
            Ckk_eu = Ckk_eu(1:so(2),1:s(2)-i*so(2));
            Ckk_eu = Ckk_eu(1:end);
            map(diag_s) = Ckk_eu;
            map(diag_i) = Ckk_eu;   
        end


end

map_ref=(xu_ref'*xu_ref)/(s(1)-1);
for i=1:s(2),
    map_ref(1:i-1,i) = map_ref(i,i);
    map_ref(i,1:i-1) = map_ref(i,i);
end
ind_0 = find(map_ref==0);
map_ref(ind_0) = 1;
nmap = map./map_ref;

if opt,
    fig_h=figure;
    map3 = [nmap nmap(:,end);nmap(end,:) nmap(end,end)];
    sur_h=surface((1:s(2)+1)'*ones(1,s(2)+1),ones(s(2)+1,1)*(1:s(2)+1),map3);
    axes_h = get(sur_h,'Parent');
    set(sur_h,'LineStyle','none');
    set(axes_h,'Box','on');
    set(axes_h,'XAxisLocation','top');
    set(axes_h,'YDir','reverse');
    paso = max(round(s(2)/10),so(2));
    paso_txt = max(5,floor(round(paso/so(2))/5)*5);
    paso = paso_txt*so(2);
    set(axes_h,'XTick',1+paso-so(2)/2:paso:s(2)+so(2)/2);
    set(axes_h,'XTickLabel',num2str((paso_txt:paso_txt:so(1))'));
    set(axes_h,'YTick',1+paso-so(2)/2:paso:s(2)+so(2)/2);
    set(axes_h,'YTickLabel',num2str((paso_txt:paso_txt:so(1))'));
    ind=[0:0.05:0.75 0.78:0.03:0.90 0.91:0.01:0.99];
    set(fig_h,'Colormap',[([ind ones(1,34)])' ([ind ones(1,4) ind(end:-1:1)])' ([ones(1,34) ind(end:-1:1)])'])

    caxis([-1 1]);

    colorbar
    tit_h=get(axes_h,'Title');
   % set(tit_h,'String',tit);
    %set(tit_h,'FontSize',18);
    axis([1,s(2)+1,1,s(2)+1]);

    cprintMV(console,'',[],-1);
end
