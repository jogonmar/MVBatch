function [Pk,Wk,Tk,Rk] = aHPCA(X,indX,pc,dk,tol)

% Adaptive Hierarchical PCA for batch process monitoring 
% Chemometr. Intell. Lab., 41:73-81, 1998.
%
% [Tb,Pb,Wt,Tt]=aHPCA(X,indX,pc,dk,tol)     % complete call
%
% INPUTS:
%
% X: (IxKJ) batch-wise unfolded 2-way batch data array (augmented data-block), I(batches) x K(sampling times)*J(variables) 
%
% indX (K x 2) = start and end variable index for each one of blocks
%
% pc (K x 1) number of principal components extracted per block
%
% dk (K x 1) weighting factor for each one of blocks
%
% tol k(1 x 1) tolerance for convergence (default 1e-09)
%
% OUTPUTS:
%
% Pk (KJ x pc) block loading matrix
%
% Wk (K x pc) super weight matrix
%
% Tk (IK x pc) super score matrix
%
% Rk (I x Kpc) block score matrix, [t1-block-1, t1-block-2,...,t1-block-K ,..., tpc-block-1, tpc-block-2, ..., tpc-block-K]
%
% codified by: Jose M. Gonzalez Martinez (jogonmar@gmail.com)
% version: 0.0
% last modification: 15/Jun/12.
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
if ndims(X)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(X);
if find(s<1), error('Incorrect content of x.'); end;
sInd = size(indX);
dmin = min(s(1),s(2)/size(indX,1));
indx = find(pc > dmin);
if indx>0, pc(indx)=dmin; end;
if nargin < 5, tol = 1e-09; end

% Initializing parameters and matrices

maxIter = 5000;
[maxPC,maxb] = max(pc);
Rk = zeros(s(1),maxPC*sInd(1));
Pk = zeros(s(2),maxPC);
Wk = zeros(2*sInd(1),maxPC);
Tk = zeros(sInd(1)*s(1),maxPC);

% Computation

for a=1:maxPC
 
    [p,t] = pca(X(:,indX(1,1):indX(1,2)),'NumComponents',1,'Centered',0);
    Tk(1:s(1),a) = t;
    Pk(indX(1,1):indX(1,2),a) = p;
  
    for k=2:sInd(1)    
        
        ik = (k-1)*s(1)+1:k*s(1);
        ikold = (k-1-1)*s(1)+1:(k-1)*s(1);
        kj = indX(k,1):indX(k,2);
        
        Tk(ik,a) = Tk(ikold,a);
        t_old = Tk(ik,a)*1000;
        iter = 0;
      
        while (sum((t_old - Tk(ik,a)).^2) > tol) & (iter < maxIter)
            Rck = [];
            iter = iter + 1;
            
            coli = (a-1)*sInd(1)+k;
            Pk(kj,a) = X(:,kj)'*Tk(ik,a);
            Rk(:,coli) = X(:,kj)*Pk(kj,a);
            Rk(:,coli) = dk*Rk(:,coli)/norm(Rk(:,coli));
            Rck = [Tk(ikold,a) Rk(:,coli)];
            Wk((k-1)*2+1:k*2,a) = Rck'*Tk(ik,a);
            t_old = Tk(ik,a);
            Tk(ik,a) = Rck*Wk((k-1)*2+1:k*2,a);
            Tk(ik,a) = Tk(ik,a)/norm(Tk(ik,a));
        end

        if iter == maxIter
            s = ['WARNING: maximum number of iterations (' num2str(maxiter) ') has been reached before convergence'];
            disp(s)
        end
    end

    % DEFLACTION
    Xaux = [];
    for k=1:sInd(1)
       ik = (k-1)*s(1)+1:k*s(1);
       kj = indX(k,1):indX(k,2);
       Xaux = [Xaux X(:,kj) - Tk(ik,a)*Pk(kj,a)'];
    end
    clear X; X = Xaux; clear Xaux;
    
end