function band = updateBound(e_pos,k,band,Kref)

% Adaptation of the upper and lower boundaries when the endpoint e* at the sampling time point k
% lies in one of the extremes of the search space. 
% The original paper is: 
% González et al. Real-time synchronization of batch trajectories for on-line multivariate statistical
% process control using Dynamic Time Warping, Chemometrics and Intelligent Laboratory Systems, 105 (2011) 195-206).
%
% CALLS:
%           band = updateBound(e_pos,k,band,Kref)   % complete call
%
%
% INPUTS:
%
% e_pos:  (1x1) k-th time point of the reference batch in which the minimum 
%           cumulative weight distance Di,j is found.
%
% k:      (1x1) current sampling time point.
%
% band:   (max(Kn)x 2) upper and lower boundaries.
%
% Kref:   number of sampling time point of the reference batch.
%
% OUTPUTS: 
%
% band:   (Kband x 2) adapted upper and lower boundaries.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)     
% last modification: Mar/2010
%
% Copyright (C) 2016  Technical University of Valencia, Valencia
% Copyright (C) 2016  José M. González Martínez
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

%% Arguments checking
routine=dbstack;
assert (nargin >= 4, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if size(band,2) ~= 2 , error('Band is not accurate for the synchronization of these multivariate trajectories.'); end

nElem = size(band,1);
incrL = 0;
incrU = 0;
decrease = false;

if k+1 < nElem
    
    if e_pos == band(k,1) && e_pos-1 >=1
        decrease = true; 
        
    else if e_pos == band(k,2) && e_pos+1<= band(nElem,2),
            incrU = 2; 
         
        end
    end

    % Boundary adaptation
    if decrease==true
         incrL = band(k+2,1) - band(k,1);      
         band(k+1:k+2,1)=band(k,1);         
    end
        
        numElem = nElem - (k+3) + 1;
        
        band(k+3:nElem,1)=max(ones(numElem,1), band((k+3):nElem,1) -incrL);
        band(k+1:nElem,2)=min(ones(numElem+2,1)*Kref, band((k+1):nElem,2)+ incrU);
       
else
    if e_pos ~= band(k,1), incrL = 1; end
      
     band = [band; [max(1,band(nElem,1)+incrL) min(band(nElem,2)+incrU,band(nElem,2))]];
end
