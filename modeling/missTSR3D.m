function rec = missTSR3D(x,pc,lag,prep,perc,iter,conv)

% Missing data imputation with Trimmed Scores Regression.
%
% rec = missTSR3D(x,pc) % impute with batch-wise unfolded model
% rec = missTSR3D(x,pc,lag,prep) % general unfolded model call
% rec = missTSR3D(x,pc,lag,prep,perc,conv,iter) % complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% pc: (1x1) number of principal components for the D-statistic.
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to 
%   the current one in the row of the unfolded matrix (0 by default).
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: trajectory centering (average trajectory subtraction)
%       2: 1 + trajectory-scaling (scales data so that each pair variable and 
%           sampling time has variance 1) (default)  
%       3: 1 + variable-scaling (scales data so that each variable has
%           variance 1)
%       4: variable centering (subtraction of the average value of each
%           variable)
%       5: 4 + variable-scaling. 
%
% perc: (1x1) maximum percentage of missing values in a row or column. (0.3
%   by default)
%
% iter: (1x1) maximum number of iterations (100 by default).
%
% conv: (1x1) convergence threshold. (1e-5 by default)
%
%
% OUTPUTS:
%
% rec: (KxJxI) recovered batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/Aug/09
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

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of pc.'); end;

if nargin < 3, lag = s(1)-1; end;
if lag<0, error('Incorrect value of lag.'); end;
if lag>s(1)-1, lag=s(1)-1; end;

if nargin < 4, prep = 2; end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;

if nargin < 5, perc = 0.3; end;
if (perc<0||perc>1), error('Incorrect value of perc.'); end;

if nargin < 6, iter = 100; end;
if (iter<1), error('Incorrect value of iter.'); end;

if nargin < 7, conv = 1e-5; end;
if (conv<0), error('Incorrect value of conv.'); end;

% Main code

s=size(x);
ind_nan=find(isnan(x));
ind_an=1-ind_nan;
xu=unfold(x,lag);
ind_nan_unf=find(isnan(xu));

nnan = unfold(isnan(x),lag);
sn = size(nnan);

% if ~isempty(find(perc*sn(1)<sum(nnan,1))),
%     error('Some columns of the unfolded matrix present too much missing values.');
% end
% if ~isempty(find(perc*sn(2)<sum(nnan,2))),
%     error('Some rows of the unfolded matrix present too much missing values.');
% end

x(ind_nan) = nan;
    
loop_ini = true;    
e0=Inf;
num0=iter;
ax = x; 
while e0>conv && num0 > 0,

    [xce,m,d] = preprocess3D_nan(ax,prep);
    
    med = zeros(s);
    for j=1:s(3),
        med(:,:,j) = m;
    end
    
    dev = zeros(s);
    for j=1:s(3),
        dev(:,:,j) = d;
    end
   
    if loop_ini,
        loop_ini = false;
        xce(ind_nan) = 0;
        x(ind_nan) = med(ind_nan);
    end
    
    ax = unfold(xce,lag);
    sa = size(ax);
    
    T=[];
    P=[];
    ax2=ax;
    for j=1:pc, % PCA model
        e1=Inf;
        num1=iter;
        t=ax(:,1);
        while e1>conv && num1 > 0,
            p=t'*ax2;
            p=p/sqrt(p*p');
            t2=ax2*p';
            e1=sum((t-t2).^2);
            t=t2;
            num1 = num1-1;
        end
        
        ax2=ax2-t*p;
        T=[T t];
        P=[P,p'];
    end
         
    if pc>0,
        T=ax*P;
        for i=1:sa(2);
            TT=ax(:,[1:i-1 i+1:end])*P([1:i-1 i+1:end],:);
            r=inv(TT'*TT)*TT'*T;
            ax2(:,i)=TT*r*P(i,:)';
        end

        rec_mod = ax;
        rec_mod(ind_nan) = ax2(ind_nan);
    else
        rec_mod = xce;
    end
    
    x2 = fold(rec_mod,s(3),lag).*dev + med;
    
    e0=sum(((x(ind_nan)-x2(ind_nan))./dev(ind_nan)).^2)/length(ind_nan);
    num0 = num0-1;
    x(ind_nan)=x2(ind_nan);
    
    ax=x;
end

rec = x;