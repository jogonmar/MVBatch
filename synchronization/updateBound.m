function band = updateBound(e_pos, k, band,Kref)

% updateBound function allows us to update the upper and lower boundaries 
% in the batch synchronization in case that the endpoint e* at the k current sampling 
% time is lying on one extreme of the band. When the duration of the ongoing batch
% is longer than the longest duration of NOC batches, with every new
% measurment coming in, the function will be extend the bands  by taking
% the value of the latter one belonging to both the upper and lower
% boundaries.
%
% INPUTS:
%
% e_pos:  k-th time point of the batch reference Bref where the minimum 
%           cumulative weight distance Di,j occurs at the current k sampling time.
%
% k:      current sampling time.
%
% band:   (max(Kn)x 2) array containing the upper and lower boundaries.
%
% Kref:   number of sampling times of the batch reference Bref.
%
% OUTPUTS: 
%
% band:   (Kband x 2) array containing the updated upper and lower
%           boundaries.
%
% CALLS:
%           band = updateBound(e_pos, k, band,Kref)   % complete call
%
% Updating function proposed by González et al. [Real-time synchronization 
% of batch trajectories for on-line multivariate statistical process
% control using Dynamic Time Warping, Chemometrics and Intelligent
% Laboratory Systems, 105 (2011) 195-206)].

% codified by: Jose Maria Gonzalez-Martinez.
% version: 1.0
% last modification: March 2010.

%% Parameters checking

if nargin < 4, error('Incorrect number of arguments.');end
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

    % Boundarries updating
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
