function [allWARPS,r] = cvBands(X,W,band,ref,warping,zeta,handle,hObject)

% Off-line cross-validation procedure to estimate the performance of the RGTW
% algorithm using different window widths. 
% The original paper is: González et al. Real-time synchronization of batch 
% trajectories for on-line multivariate statistical process control using 
% Dynamic Time Warping, Chemometrics and Intelligent % Laboratory Systems, 
% 105 (2011) 195-206).
%
% CALLS:
%
%       [allWARPS,r] = cvBands(X,W,band,ref,warping)                      % no output in GUI object and default zeta
%       [allWARPS,r] = cvBands(X,W,band,ref,warping,zeta)                 % no output in GUI object 
%       [allWARPS,r] = cvBands(X,W,band,ref,warping,zeta, handle,hObject) % complete call
%
%
% INPUTS:
%
% X: (1xI) cell array containing the measurements collected for J variables at 
%    Ki different sampling times for each one of the I batches. 
%       
% W: (JxJ) matrix containing weights to give more importance to certain
%    variables based on a criterium selected in the offline synchronization.
%
% band: (max(Ki)x 2) matrix containing the upper and lower limits that define
%       the research space to estimate the local and cumulative distances, and
%       hence, the warping path. Note that max(Ki) is the maximum batch duration of
%       I historical batches.
%
% ref: index of the batch selected as reference for batch synchronization.
%
% warping: (1xI) cell array containing the warping information from the
%           offline synchronization of I historical batches.
%
% zeta: number of windows widths to be used to study the performance of the
%       RGTW.
%
% handle: GUI handle to access the text box to report process on CV.
%
% hObject: GUI object to update progress on the CV synchronization.
%
%
% OUTPUTS:
%
% allWARPS: (1xI) cell array containing the warping information obtained
%           for each one of the stated number of window width using the RGTW algorithm.
%
% r:        (zeta x I) Pearson's correlation factors calculated from the
%           warping information using DTW and RGTW, where zeta is the
%           number of window widths taken into account in the study and I
%           the number of batches synchronized.
%
%
% coded by: José M. González Martínez (J.Gonzalez-Martinez@shell.com)
% last modification: May/10
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

%% Parameters checking and settings

if nargin < 5, error('Incorrect number of input paramters. Please, check the help for further details.'); end
if ~iscell(X), error('The data structure must be a cell array to contain the unsynchronized batch trajectories.'); end
sw = size(W);
if sw(1) ~= sw(2), error('The weight matrix must be a squared one.'); end
if min(W) < 0, error('Matrix W must be positive definite'); end
nJ = size(X{1},2);
s =size(X);
for i=1:s(2)
    if size(X{i},2)~=nJ, error('Not all batches have the same number of variables collected.'); end
end
if nJ~=sw(1), error('The number of the weights do not match with the number of process variables.'); end
if size(band,2)~= 2, error('Dimensions of the array not expected. '); end
if ref < 1 || ref > s(2), error('Wrong selection of the reference batch.'); end
if length(warping) ~= s(2), error('The number of cells containing the warping information must be the same as the number of batches.'); end
if nargin < 6, zeta = 4; end
if zeta < -1, error('The number of window width to be studied must be greater than 0'); end
plotting = true;
if nargin < 7, plotting = false; end
if nargin == 7 && nargin < 8, error('The object handle is required to update the synchronization information in the user interface.'); end

% Parameter initialization
allWARPS = cell(1,zeta);
r = zeros(zeta,size(X,2));

if plotting
    set(handle.uite_DTW_Window,'String','');
    guidata(hObject, handle);
end
pause(.1);

%% Iterative cross-validation procedure

for w=1:zeta
 
    if plotting
        cprint(handle.uite_DTW_Window, strcat(strcat('CV with window width = ',num2str(w))));
        pause(.01);
        guidata(hObject, handle);
    end

    warpOri = cell(1,size(X,2));
    WARP = zeros(size(X{ref},1),s(2));
    
    for i=1:s(2)
        % Split up the X data matrix in a calibration and test set
        numItems=1;
        train=  cell(1,s(2)-1); 

        test = X{i};

        if(i-1)>=1
           for j=1:i-1 
               train{numItems}=X{j};
               numItems = numItems + 1;
           end        
        end

        if(i+1)<=s(2)
            for j=i+1:s(2)
                train{numItems}=X{j};
                numItems = numItems + 1;
            end
        end

        % On-line synchronization using a window width w. 

        %Inner validation

        [trainSc,rngTr] = scale_(train);

        Bref = scale_(X{ref},rngTr);

        %warp = cell(1,s(2)-1);
        warp = zeros(size(X{ref},1),s(2)-1);
        for b=1:s(2)-1
            [sBn,warp(:,b)] = onSyn(train{b},Bref, band,W,w,rngTr);
        end

        maxi = 0; mini = 9e10;
        for j=1:s(2)-1
            maxi = max(maxi,size(train{j},1));
            mini = min(mini,size(train{j},1)); 
        end

        % Estimate the synchronization band
        bandInner = estimationBD(warp);

        % Outer validation
        [sBn,WARP(:,i) warpOnSyn] = onSyn(test,Bref,bandInner,W,w,rngTr);

        % PEARSON'S CORRELATION FACTOR ESTIMATION %

         % Calculate the Pearson's correlation factor between the off-line and
         % on-line warping information estimated for a window width w.

         items = warpOnSyn(size(warpOnSyn,1),2); 
         indxVal = zeros(1,items);
         indxOri = zeros(1,items);

         for k=1:items
             indxOri(k) = warping{i}(max(find(warping{i}(:,2)==k)),1);
             indxVal(k) = warpOnSyn(max(find(warpOnSyn(:,2)==k)),1);
         end    
         r(w,i)= corr(indxOri',indxVal');

         if plotting
             cprint(handle.uite_DTW_Window, strcat(strcat('Leaving batch #',num2str(i)),' out'));
             guidata(hObject, handle);
         end
    end

allWARPS{w} = WARP;

end

%--------------------------------------------------------------------------

