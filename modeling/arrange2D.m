function batchesO = arrange2D(batchesI,inter,units,method),

% Normalizes the sampling interval inter and intra variables in a number of
% batches.
%
% batchesO = arrange2D(batchesI)    % normalize to greatest common divisor
% batchesO = arrange2D(batchesI,inter,units,method)    % complete call
%
% INPUTS:
%
% batchesI: (1xI) group of batches with a common set of variables (some 
%   of the variables may not have been measured in some batches).
%   
%   batchesI(i).data: cell (1xM) data corresponding to the M different sampling
%       timing in batch i. 
%
%       batchesI(i).data{m}: (Kx(1+J)) data corresponding to each different 
%           sampling timing (m) in a batch (i). The matrix contains in the 
%           first column the sampling timing and in the rest the values 
%           collected for J variables according to that sampling timing.
%
% inter: (1x1) normalization method.
%       0: Obtain the most common sampling interval of each variable and
%           normalize to their greatest common divisor. (default)
%       1: Obtain the most common sampling interval of each variable and
%           normalize to their minimum.
%       2: Obtain the most common sampling interval of each variable and
%           normalize to their maximum. 
%       3: Obtain the most common sampling interval of each variable and
%           normalize to their least common multiple.
%
% units: (1x1) units used for sampling interval rounding in methods 0 and 3
%   (1e-3 by default).
%
% method: see methods in function interp1:
%       'nearest': nearest neighbor interpolation (default)
%       'linear': linear interpolation
%       'spline': piecewise cubic spline interpolation (SPLINE)
%       'pchip': shape-preserving piecewise cubic interpolation
%       'cubic': same as 'pchip'
%       'v5cubic': the cubic interpolation from MATLAB 5, which does not
%                   extrapolate and uses 'spline' if X is not equally
%                   spaced.
%
%
% OUTPUTS:
%
% batchesO: cell (1xI) group of batches with a common set of variables (some 
%   of the variables may not have been measured in some batches).
%
%   batchesO{i}: (Kix(1+J2)) data corresponding to the original time ordering and the 
%       J2 variables collected in the batch with a common sampling timing. 
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
s = size(batchesI);
if s(1)~=1, error('Incorrect content of batchesI.'); end; 
for i=1:s(2),
    if ~isfield(batchesI(i),'data'),
        error('Incorrect content of batchesI.'); 
    elseif ~iscell(batchesI(i).data),
        error('Incorrect content of batchesI.'); 
    else
        s2 = size(batchesI(i).data);
        if s2(1)~=1, error('Incorrect content of batchesI.'); end; 
        for j=1:s2(2),
            s3 = size(batchesI(i).data{j});
            if s3(2)==1, error('Incorrect content of batchesI.'); end; 
        end
    end;
end

if nargin < 2, inter = 1; end;
if (inter<0||inter>3), error('Incorrect value of inter.'); end;
if nargin < 3, units = 1e-3; end; 
if nargin < 4, method = 'spline'; end; 

% Delete non-measured variables

for i=1:length(batchesI), 
    for o=1:length(batchesI(i).data), 
        if ~isempty(batchesI(i).data{o})
           var{o}.lot{i}=batchesI(i).data{o}(:,1);
        end
    end
end

ex_var=ones(1,length(var));
for o=1:length(var),
    if isempty(var{o}),
        ex_var(o) = 0;
    end
end;

for i=1:length(batchesI),
    batchesI(i).data = batchesI(i).data(find(ex_var));
end            
var = var(find(ex_var));  

% Compute the time interval for each variable

for o=1:length(var),
    inter2=[];
    for i=1:length(var{o}.lot),
        if ~isempty(var{o}.lot{i}),
           int = var{o}.lot{i}(2:end) - var{o}.lot{i}(1:end-1); 
           inter2 =[inter2;int];
        end
    end
    rate(o) = mode(inter2); 
end


% Interpole and compute the measurement vector

switch inter,
    case 0, % Greatest common divisor
        rate = round(rate/units);
        rate_tot=rate(1);
        for o=2:length(rate)
            rate_tot=gcd(rate_tot,rate(o));
        end
        
    case 1, % Interpole to maximum rate
        rate = rate/units;
        rate_tot=min(rate);
           
    case 2, % Interpole to minimum rate
        rate = rate/units;
        rate_tot=max(rate);
    
    case 3, % Least common multiple
        rate = round(rate/units);
        rate_tot=rate(1);
        for o=2:length(rate)
            rate_tot=lcm(rate_tot,rate(o));
        end
        
    otherwise,
        exit('Unknown normalization method.')
        
end

for i=1:length(batchesI), 
   ini=-Inf;
   fin=Inf;
   batchesO{i}=[];
   
   for o=1:length(batchesI(i).data),
       if ~isempty(batchesI(i).data{o}),
           if batchesI(i).data{o}(1,1)>ini,
               ini = batchesI(i).data{o}(1,1);
           end
                    
           if batchesI(i).data{o}(end,1)<fin,
               fin = batchesI(i).data{o}(end,1);
           end
       end
   end

   first = true;
   for o=1:length(var),
       if ~isempty(batchesI(i).data{o}),                    
           ind_i = find(batchesI(i).data{o}(:,1)>=ini);           
           ind_f = find(batchesI(i).data{o}(:,1)>=fin);
                    
           orig = batchesI(i).data{o}(ind_i(1):ind_f(1),1:end);
           data = interp1(orig(:,1),orig,(ini:rate_tot*units:fin)',method,'extrap');
           
           if first,         
                batchesO{i} = [batchesO{i} data];
                first = false;
           else
                batchesO{i} = [batchesO{i} data(:,2:end)];
           end
       else
           batchesO{i} = [batchesO{i} nan*ones(length(ini:rate_tot*units:fin),1)];
       end
   end
end