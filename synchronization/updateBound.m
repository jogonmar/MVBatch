function band = updateBound(e_pos,k,band,Kref)

% Update upper and lower boundaries in real-time batch synchronization.
% The original work is:
% [1] González et al. Real-time synchronization of batch trajectories for on-line 
% multivariate statistical process control using Dynamic Time Warping, Chemometrics 
% and Intelligent Laboratory Systems, 105 (2011) 195-206).
% [2] González-Martinez, J.M. Advances on bilinear modeling of biochemical
% batch processes (2015). PhD thesis, DOI: 10.4995/Thesis/10251/55684.
%
% CALLS:
%           band = updateBound(e_pos,k,band,Kref)   % complete call
%
% INPUTS:
%
% e_pos:  k-th time point of the reference batch where the minimum 
%           cumulative weight distance Di,j occurs at the current k sampling time.
%
% k:      current sampling time.
%
% band:   (max(Kn)x 2) array containing the upper and lower boundaries.
%
% Kref:   number of sampling times of the reference batch.
%
% OUTPUTS: 
%
% band:   (Kband x 2) updated upper and lower boundaries.
%
% codified by: Jose Maria Gonzalez-Martinez.
% last modification: Mar/10.
%
% Copyright (C) 2016  José M. Gonzalez-Martinez
% Copyright (C) 2016  Technical University of Valencia, Valencia
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

%% Parameters checking

if nargin < 4, error('Incorrect number of input paramters. Please, check the help for further details.');end
if size(band,2) ~= 2 , error('Dimensions of the array not expected.'); end

%% Computation
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

    % Bound adaptation
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
