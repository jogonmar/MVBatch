function [ncpr,npr,cpr,pr] = mv_parreg(x,lag,type,p_min,opt,cumul,console)

% Computes the partial regression function of a measurement vector..
%
% [cpr,pr] = mv_parreg(x)          % standard inputs
% [cpr,pr] = mv_parreg(x,lag,type,p_min,opt,cumul) % output in MATLAB console 
% [cpr,pr] = mv_parreg(x,lag,type,p_min,opt,cumul,console)  % complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
% one in the row of the unfolded matrix (10 by default), minimum 1 and maximum K-1.
%
% type: (1x1) type of covariance map.
%       1: dynamic partial covariance map (pseudoinverse) (by default). 
%       2: instantaneous+dynamic partial covariance map (pseudoinverse). 
%
% p_min: (1x1) value in the interval (0,1] to set the minimum portion of
%   the total variance for an eigenvalue to be taken into account (0.1 by default).
%
% opt: boolean (1x1) 
%       true: plot results. (by default)
%       false: do not plot results.
%
% cumul: boolean (1x1) 
%       true: plot the cumulative. (by default)
%       false: plot the variables sepparately.  
% 
%
% console: (1x1) handle of the EditText of the interface, 0 stands for the
%   MATLAB console (by default)
%
%
% OUTPUTS:
%
% ncpr: (1x(lag+1)) normalized cummulative partial regression sum of squares per lag.
%
% npr: (Jx(lag+1)) normalized partial regression sum of squares per variable and lag.
%
% cpr: (1x(lag+1)) cummulative partial regression sum of squares per lag.
%
% pr: (Jx(lag+1)) partial regression sum of squares per variable and lag.
%
%
% codified by: José Camacho Páez.
% last modification: 08/Sep/16.

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 2, lag = 10; end;
if lag<1, error('Incorrect value of lag.'); end;
if lag>s(1)-1, lag=s(1)-1; end;
if nargin < 3, type = 1; end;
if (type<1||type>2), error('Incorrect value of type.'); end;
if nargin < 4, p_min = 0.1; end;
if (p_min<0||p_min>1), error('Incorrect value of p_min.'); end;
if nargin < 5, opt = 1; end;
if nargin < 6, cumul = 1; end;
if nargin < 7, console = 0; end;

% Computation

xu = unfold(x,lag);
s = size(xu);
so = size(x);
  
pr = zeros(so(2),lag);
cpr = zeros(1,lag);

cprintMV(console,'Processing.... Please, be patient.',[],0);
               
switch type,
 
    case 1, 
        eu = xu; % backward regression
        au = xu; % fordward regression
            
        for i=1:lag, % time distance in MV units
            % Backward regression
            for o=s(2):-so(2):(i+1)*so(2),
                indk = o-(i+1)*so(2)+1:o-i*so(2);
                indk1 = o-so(2)+1:o;
                p = au(:,indk)'*au(:,indk); % Pseudo-inverse with specific tolerance
                est_eu(:,indk1) = au(:,indk)*(pinv(p,p_min*trace(p))*(au(:,indk)'*eu(:,indk1)));                 
            end
            
            Ckk_eu = (est_eu(:,s(2)-so(2)+1:s(2))'*est_eu(:,s(2)-so(2)+1:s(2)))/(s(1)-1);
            pr(:,i) = diag(Ckk_eu);
            cpr(i) = sum(pr(:,i))/so(2);
            
            % Forward regression
            for o=1:so(2):s(2)-(i+1)*so(2),
                indk = o+i*so(2):o+so(2)*(1+i)-1;
                indk1 = o:o+so(2)-1;
                p = eu(:,indk)'*eu(:,indk); % Pseudo-inverse with specific tolerance
                au(:,indk1) = au(:,indk1) - eu(:,indk)*(pinv(p,p_min*trace(p))*(eu(:,indk)'*au(:,indk1)));
            end
            
            eu(:,(i+1)*so(2):s(2)) = eu(:,(i+1)*so(2):s(2)) - est_eu(:,(i+1)*so(2):s(2));
        end

    case 2, 
        eu = xu; % backward regression
        au = xu; % fordward regression
        
        o=s(2)-so(2)+1;
        for i=0:so(2)-1, % instantaneous regression
            ind_rest = [o:o+i-1 o+i+1:o+so(2)-1];
            ind = o+i;
            p = au(:,ind_rest)'*au(:,ind_rest); % Pseudo-inverse with specific tolerance
            est_eu(:,i+1) = au(:,ind_rest)*(pinv(p,p_min*trace(p))*(au(:,ind_rest)'*eu(:,ind)));
            ini_eu(:,i+1) = eu(:,ind) - est_eu(:,i+1);
        end   
        
        Ckk_eu  = (est_eu'*est_eu)/(s(1)-1);
        pr(:,1) = diag(Ckk_eu );
        cpr(1) = sum(pr(:,1))/so(2);
            
        for i=1:lag, % time distance 

            % Forward regression without specific variable
            indk = s(2)-so(2)+1:s(2);
            indk1 = s(2)-(i+1)*so(2)+1:s(2)-i*so(2);
                
            for j=0:so(2)-1, 
                ind_rest = indk([1:j j+2:end]);
                p = eu(:,ind_rest)'*eu(:,ind_rest); % Pseudo-inverse with specific tolerance
                est_au(:,indk1) = eu(:,ind_rest)*(pinv(p,p_min*trace(p))*(eu(:,ind_rest)'*au(:,indk1)));
                ausp(:,:,j+1) = au(:,indk1) - est_au(:,indk1);
            end
            
            for o=1:so(2):s(2)-i*so(2),
                indk = o+i*so(2):o+so(2)*(1+i)-1;
                indk1 = o:o+so(2)-1;
                p = eu(:,indk)'*eu(:,indk); % Pseudo-inverse with specific tolerance
                est_au(:,indk1) = eu(:,indk)*(pinv(p,p_min*trace(p))*(eu(:,indk)'*au(:,indk1)));
            end          

            % Backward regression without specific variable
            for o=s(2):-so(2):(i+1)*so(2),
                indk = o-(i+1)*so(2)+1:o-i*so(2);
                indk1 = o-so(2)+1:o;
  
                p = au(:,indk)'*au(:,indk); % Pseudo-inverse with specific tolerance
                est_eu(:,indk1) = au(:,indk)*(pinv(p,p_min*trace(p))*(au(:,indk)'*eu(:,indk1)));
                eu(:,indk1) = eu(:,indk1) - est_eu(:,indk1); 
            end
            
            est_eu=[];
            for j=0:so(2)-1, 
                p = ausp(:,:,j+1)'*ausp(:,:,j+1); % Pseudo-inverse with specific tolerance
                est_eu(:,1+j) = ausp(:,:,j+1)*(pinv(p,p_min*trace(p))*(ausp(:,:,j+1)'*ini_eu(:,1+j)));
                pr(1+j,i+1) = est_eu(:,1+j)'*est_eu(:,1+j)/(s(1)-1); 
            end
            ini_eu = ini_eu - est_eu; 
                
            au(:,1:s(2)-i*so(2)) = au(:,1:s(2)-i*so(2)) - est_au(:,1:s(2)-i*so(2));
 
            cpr(i+1) = sum(pr(:,i+1))/so(2);  
            
        end

end
 
tot=diag(xu(:,s(2)-so(2)+1:s(2))'*xu(:,s(2)-so(2)+1:s(2)))/(s(1)-1);
npr = pr./(tot*ones(1,lag+type-1));
tot=sum(tot)/so(2);
ncpr = cpr/tot;

cprintMV(console,'',[],-1);

if opt,
    fig_h = figure;
    if cumul,
        
        if type==2,
            b = bar(0:length(ncpr),[ncpr 1-sum(ncpr)]);
            axes_h = get(b,'Parent');
            set(axes_h,'XTickLabel',strvcat(num2str((0:lag)'),'Res')); 
        else
            b = bar(1:length(ncpr)+1,[ncpr 1-sum(ncpr)]);
            axes_h = get(b,'Parent'); 
            set(axes_h,'XTickLabel',strvcat(num2str((1:lag)'),'Res')); 
        end
        set(axes_h,'FontSize',12);
        xl_h=get(axes_h,'XLabel');
        set(xl_h,'String','LMVs');
        set(xl_h,'FontSize',14);
        yl_h=get(axes_h,'YLabel');
        set(yl_h,'String','Normalized Mean Sum of Squares');
        set(yl_h,'FontSize',14);
        tit_h=get(axes_h,'Title');    
%         set(tit_h,'FontSize',18);
        axis([0,length(cpr)+1-(type~=1),0,1]);
    else  
        s = size(pr);
        if type==2,
            b = bar3([npr 1-sum(npr,2)]);
            axes_h = get(b(1),'Parent');
            set(axes_h,'XTickLabel',strvcat(num2str((0:lag)'),'Res')); 
        else
            b = bar3([npr 1-sum(npr,2)]);
            axes_h = get(b(1),'Parent'); 
            set(axes_h,'XTickLabel',strvcat(num2str((1:lag)'),'Res')); 
        end

        set(axes_h,'FontSize',12);
        tit_h=get(axes_h,'Title');
%         set(tit_h,'FontSize',20);
        xl_h=get(axes_h,'XLabel');
        set(xl_h,'String','LMVs');
        set(xl_h,'FontSize',14);
        yl_h=get(axes_h,'YLabel');
        set(yl_h,'String','Variables');
        set(yl_h,'FontSize',14);
        zl_h=get(axes_h,'ZLabel');
        set(zl_h,'String','Normalized Sum of Squares');
        set(zl_h,'FontSize',14);     
        axis([0,length(cpr)+1.5,0.5,so(2)+0.5,0,1]);
        a=get(axes_h,'View');
        set(axes_h,'View',[-a(1),a(2)]);
    end

    if type==2,
        set(tit_h,'String','Partial Regression Sum of Squares');
    else
        set(tit_h,'String','Dynamic Partial Regression Sum of Squares');
    end
end
