function [xp,average,scale] = preprocess3D(x,prep)

% Preprocess 3-way batch data.
%
% [xp,average,scale] = preprocess3D(x)          % for trajectory scaling
% [xp,average,scale] = preprocess3D(x,prep)     % complete call
%
%
% INPUTS:
%
% x: (KxJxI) three-way batch data matrix, K(sampling times) x J(variables)
%   x I(batches)
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
%
% OUTPUTS:
%
% xp: (KxJxI) preprocessed data.
%
% average: (K x J) sample average according to the preprocessing method.
%
% scale: (K x J) sample scale according to the preprocessing method.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 14/Oct/08
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
if ndims(x)~=3, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 2, prep = 2; end;
if (prep<0||prep>5), error('Incorrect value of prep.'); end;

% Computation

% if prep>2,
%     disp('Not all the functionality is avaiable in this version')
%     prep = 2
% end;

K=1;
J=2;
I=3;

xp=zeros(s);

switch prep,
    
    case 1, % Trajectory centering
        average = sum(x,3)/s(I);    
        scale = ones(s(K),s(J));
        for j=1:s(J),
            xp(:,j,:)=squeeze(x(:,j,:))-average(:,j)*ones(1,s(I));
        end;
        
    case 2, % Trajectory centering and scaling
        xc=zeros(s);
        average = sum(x,3)/s(I);
        for j=1:s(J),
            x_red=squeeze(x(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xc(:,j,:)=x_red-average(:,j)*ones(1,s(I));
        end
        scale = sqrt(sum(xc.^2,3)/(s(I)-1));
        scale(find(scale==0))=1; 
        for j=1:s(J),
            x_red=squeeze(xc(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xp(:,j,:)=x_red./(scale(:,j)*ones(1,s(I)));
        end;
        
    case 3, % Trajectory centering and variable scaling
        xc=zeros(s);
        average = sum(x,3)/s(I);
        for j=1:s(J),
            x_red=squeeze(x(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xc(:,j,:)=x_red-average(:,j)*ones(1,s(I));
        end
        scale = ones(s(K),1)*sqrt(sum(sum(xc.^2,3),1)/((s(I)-1)*s(K)));
        scale(find(scale==0))=1; 
        for j=1:s(J),
            x_red=squeeze(xc(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xp(:,j,:)=x_red./(scale(:,j)*ones(1,s(I)));
        end;
        
    case 4, % Variable centering 
        average = ones(s(K),1)*sum(sum(x,3),1)/(s(I)*s(K));    
        scale = ones(s(K),s(J)); 
        for j=1:s(J),
            x_red=squeeze(x(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xp(:,j,:)=x_red-average(:,j)*ones(1,s(I));
        end;

    case 5, % Variable centering and variable scaling        
        xc=zeros(s);
        average = ones(s(K),1)*sum(sum(x,3),1)/(s(I)*s(K));
        for j=1:s(J),
            x_red=squeeze(x(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xc(:,j,:)=x_red-average(:,j)*ones(1,s(I));
        end
        scale = ones(s(K),1)*sqrt(sum(sum(xc.^2,3),1)/(s(I)*s(K)-1));
        scale(find(scale==0))=1; 
        for j=1:s(J),
            x_red=squeeze(xc(:,j,:));
            if size(x_red,2)==1, x_red=x_red'; end;
            xp(:,j,:)=x_red./(scale(:,j)*ones(1,s(I)));
        end;
        
    otherwise, % No preprocessing 
        average = zeros(s(K),s(J));     
        scale = ones(s(K),s(J)); 
        xp = x;
end


