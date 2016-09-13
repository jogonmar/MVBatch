function [cumpress,press,pem] = crossval3D_s(x,pc,lag,clu,leave_m,blocks_r,blocks_c,fold_m,prep,order)

% Cross-validation for square-prediction-errors computing in 3-to-2 way 
%   batch data, divisions in the samples.  
%
% [cumpress,press,pem] = crossval3D_s(x,pc) % Variable-wise data
% [cumpress,press,pem] = crossval3D_s(x,pc,Inf) % Batch-wise data
% [cumpress,press,pem] = crossval3D_s(x,pc,lag,clu) % Multi-phase data (division in the samples)                                                                                            
% [cumpress,press,pem] = crossval3D_s(x,pc,lag,clu,leave_m,blocks_r,blocks_c,fold_m,prep,order) % Complete call
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
%
% pc: (1x1) number of principal components.
%
% lag: (1x1) number of immediate lagged measurement-vectors (LMVs) added to the current
%   one in the row of the unfolded matrix (0 by default).
%
% clu: (1xK) vector with the assignment of the sampling times to the phases, numbered
%   from 1 onwards (1 phase by default). 
%
% leave_m: (text) cross-validation procedure:
%   'rkf': row-wise k-fold cross-validation.
%   'skf': sample-wise k-fold cross-validation (by default).
%   'iskf': iterative sample-wise k-fold cross-validation.
%   'cskf': cross-corrected sample-wise k-fold cross-validation. 
%
% blocks_r: (1x1) maximum number of blocks of samples (min(I,30) by default)
%
% blocks_c: (1x1) maximum number of blocks of variables (Inf by default).
%
% fold_m (text) folding method:
%   'mean': mean of all the values of a variable.
%   'first': first value (by default).
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
% order: (structure) to define a constant random ordering of columns and
%   rows.
%       order.input: (boolean)
%           true: take ordering from the structure.
%           false: compute it ramdomly (by default).
%       order.cols: (1xn_cols) columns ordering.
%       order.rows: (1xn_rows) rows ordering.
%
%
% OUTPUTS:
%
% cumpress: (1x1) cumulative PRESS.
%
% press: (Kx1) PRESS throughout the batch time.
%
% pem: (KxJxI) Matrix containing prediction errors.
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

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, lag = 0; end;

if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;

if nargin < 4, clu = ones(s(1),1); end;
if nargin < 5, leave_m = 'skf'; end;
if nargin < 6, blocks_r = min(s(3),30); end;
if nargin < 7, blocks_c = Inf; end;
if nargin < 8, fold_m = 'first'; end;
if nargin < 9, prep = 2; end;
if nargin < 10, order.input = false; end;

if ~isequal(leave_m,'skf'),
    error(-1)
    disp('Not all the functionality is avaiable in this version')
    leave_m = 'skf'
end;

if pc<0, error('Incorrect value of pc.'); end;
if lag<0, error('Incorrect value of lag.'); end;
if lag>s(1)-1, lag=s(1)-1; end;
if ndims(clu)~=2, error('Incorrect number of dimensions of clu.'); end;
if find(size(clu)~=[s(1) 1]), error('Incorrect content of clu.'); end;
if lag>1 & clu(lag)>1, error('First phase shorter than lag.'); end;
if blocks_r>s(3), blocks_r = s(3); end
if (blocks_r<2), error('Incorrect value of blocks_r.'); end;
if blocks_c>s(2)*(lag+1), blocks_c = s(2)*(lag+1); end
if (blocks_c<2), error('Incorrect value of blocks_c.'); end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;


% Initialization

pem = zeros(size(x));
press = zeros(s(1),1);
cumpress = 0;

if ~order.input,
    rows = rand(1,s(3));
else
    rows = order.rows;
end

[a,r_ind]=sort(rows);
elem_r=s(3)/blocks_r;

% Cross-validation
for i=1:blocks_r,
    ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection %COMPROBAR
    i2 = ones(s(3),1);
    i2(ind_i)=0;
    sample = x(:,:,ind_i);
    calibr = x(:,:,find(i2)); 
    
    [ccs,av,st] = preprocess3D(calibr,prep);
    
    sc=sample;
    scs=sample;
    for j=1:length(ind_i),
        scs(:,:,j) = (sc(:,:,j)-av)./st;
    end
            
    if pc > 0, % PCA Modelling
        m = max(clu);
        srec_3D = zeros(size(sample));
            
        for o=1:m, % Phases
            indx = find(clu==o);
            indx2 = [(max(1,indx(1)-lag):indx(1)-1)';indx];
         
            c_2D =unfold(ccs(indx2,:,:),lag);
            sizc=size(c_2D);
            s_2D = unfold(scs(indx2,:,:),lag);
            sizs=size(s_2D);
            
            avc_prep=ones(sizc(1),1)*mean(c_2D);
            avs_prep=ones(sizs(1),1)*mean(c_2D);
            
            if ~order.input,
                cols = rand(1,sizs(2));
            else
                cols = order.cols(1:sizs(2));
            end            
            [a,c_ind]=sort(cols);
            elem_c=sizs(2)/blocks_c;
            
            srec_2D=[];
            pem_2D=[];    
            switch lower(leave_m)
    
                % Leave-n-batches-out cross_validation
                case 'rkf',
                    p =  pcamv(c_2D,'Centered',false,'NumComponents',pc);
                    t_est = s_2D*p;
                    srec_2D = t_est*p';                   
                    srec_3D = scs(indx2,:,:)-fold(srec_2D,length(ind_i),lag,fold_m); 
                    
                % Leave-n-samples-out cross_validation based on zero
                % values.
                case 'skf',
                    p =  pcamv(c_2D,'Centered',false,'NumComponents',pc);
                    t_est = s_2D*p;
                    srec_2D = t_est*p';
                    erec_2D = s_2D - srec_2D;
                    for j=1:blocks_c,
                        ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection                                  
                        pem_2D(1:sizs(1),ind_j) = (s_2D(:,ind_j)-avs_prep(:,ind_j))*(p(ind_j,:)*p(ind_j,:)') + erec_2D(:,ind_j);
                    end
                    srec_3D = fold(pem_2D,length(ind_i),lag,fold_m);

                % Leave-n-samples-out cross_validation based on zero
                % values until convergence.
                case 'iskf',
                    p =  pcamv(c_2D,'Centered',false,'NumComponents',pc);
                    t_est = s_2D*p;
                    srec_2D = t_est*p';
                    erec_2D = s_2D - srec_2D;
                    for j=1:blocks_c,
                        ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection    
                        pem_2D(1:sizs(1),ind_j) = (inv(eye(length(ind_j))-p(ind_j,:)*p(ind_j,:)') * erec_2D(:,ind_j)')';
                    end
                    srec_3D = fold(pem_2D,length(ind_i),lag,fold_m);
                
                % Cross-valiadtion corrected-leave-n-samples-out cross_validation based on zero
                % values 
                case 'cskf',

                    [p,t] =  pcamv(c_2D,pc);
                    [rec,av] = preprocess2D(t);

                    rec_sam=s_2D*p;
                    for j=1:length(ind_i),
                        rec_sam(j,:) = (rec_sam(j,:)-av);
                    end

                    p2 =  pcamv([c_2D rec],pc);
                    s_2D2 = [s_2D rec_sam];
                    t_est = s_2D2*p2;
                    srec_2D = t_est*p2';
                    erec_2D = s_2D2 - srec_2D;
                    for j=1:blocks_c,                    
                        ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection
                        pem_2D(1:sizs(1),ind_j) = (s_2D2(:,ind_j)-avs_prep(:,ind_j))*(p2(ind_j,:)*p2(ind_j,:)') + erec_2D(:,ind_j);
                    end
                    
                    srec_3D = fold(pem_2D,length(ind_i),lag,fold_m);
                    
                otherwise
                    error('Incorrect leave_m.');

            end
              
            ini=find(indx2==indx(1));
            pem(indx,:,ind_i) = srec_3D(ini:end,:,:);
                 
            if isequal(lower(fold_m),'mean'), % Error correction between phases
                if indx(1)>lag+1,
                    w_ind=(indx(1)-lag:indx(1)-1);    
                    w1=(length(w_ind):-1:1)'*ones(1,s(2));
                    w2=(1:length(w_ind))'*ones(1,s(2));
                    w2=min(w2,length(indx));
                    
                    for j=1:length(ind_i); 
                        pem(w_ind,:,ind_i(j)) = (w1.*pem(w_ind,:,ind_i(j)) + w2.*srec_3D(1:ini-1,:,j))./(w1+w2);
                    end
                end
            end
        end 

    else % Modelling with the average
        pem(:,:,ind_i) = scs;
    end

end

press = sum(sum(pem.^2,3),2);

cumpress = sum(press(find(clu)));
