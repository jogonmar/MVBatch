function [Sout Rout] = S_pairwise(Xin)
% Missing Data Imputation Toolbox v1.0
% A. Folch-Fortuny, F. Arteaga and A. Ferrer
% Copyright (C) 2015 F. Arteaga
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
%
%
% INPUTS
%
% Xin: data matrix with missing data
%
%
% OUTPUTS
%
% Sout: covariance matrix
% Rout: correlation matrix




[n p]=size(Xin);

nonan=ones(n,p)-isnan(Xin);

% Variances estimation

% h = waitbar(X,'message', property, value, property, value, ...) creates and displays a waitbar of fractional length X. The handle to the waitbar figure is returned in H. X should be between 0 and 1.
h = waitbar(0,'Please wait...',...
    'Color',[0.5 0.6 0.7],... % background color
    'Name','Variances estimation...',... % figure title
    'Position',[500 300 364 60]); % position of the figure
waitbar(0/p,h,[num2str(ceil(0/p*100)) '%'])

for i=1:p,
  x=[];
  for j=1:n,
    if nonan(j,i)==1,
      x=[x;Xin(j,i)];
    end
  end
  Sout(i,i)=var(x);
  

    waitbar(i/p,h,[num2str(ceil(i/p*100)) '%'])

  
end

close(h)



% Covariance estimation

% h = waitbar(X,'message', property, value, property, value, ...) creates and displays a waitbar of fractional length X. The handle to the waitbar figure is returned in H. X should be between 0 and 1.
h = waitbar(0,'Please wait...',...
    'Color',[0.5 0.6 0.7],... % background color
    'Name','Covariances estimation...',... % figure title
    'Position',[500 300 364 60]); % position of the figure
waitbar(0/(p-1),h,[num2str(ceil(0/p*100)) '%'])

for i=1:p-1,
  for j=i+1:p, % for each pair Xi,Xj
    xy=[];
    for k=1:n, % considers only the rows with data in both Xi and Xj
      if nonan(k,i)+nonan(k,j)==2,
        xy=[xy;Xin(k,i) Xin(k,j)];
      end
    end
    vxy=cov(xy);
    if numel(vxy)>1
        Sout(i,j)=vxy(1,2);
    else
        Sout(i,j)=0;
    end
    Sout(j,i)=Sout(i,j);
  end


  
  S_aux=Sout;
  K=size(S_aux,1);
  vari=ones(1,K);

  for k=1:K-1
    for j=i+1:K
      Rout(k,j)=S_aux(k,j)*sqrt(vari(k)*vari(j)/(S_aux(k,k)*S_aux(j,j))); Rout(j,k)=Rout(k,j);
    end
  end

  for k=1:K
    Rout(k,k)=vari(k); 
  end
  

    waitbar(i/(p-1),h,[num2str(ceil(i/p*100)) '%'])

    
end

close(h)