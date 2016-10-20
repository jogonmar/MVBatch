function Bref = refBatch(X,method)

% Find a reference batch based on a specific criterium.
%
% CALLS: 
%           refBatch(X,method)      % complete call
%
%
% INPUTS:
%
% X: (1xI) cell array containing the measurements collected for J variables at 
%    Ki different sampling times for each of the I batches. 
%           
% method:   criterium to be followed to estimate the batch reference (median or average).
%           
% OUTPOUT:
%
% Bref:     batch ID selected to play the role of the reference batch in batch synchronization.
%
% coded by: José M. Gonzalez-Martinez (J.Gonzalez-Martinez@shell.com)          
% last modification: May/11.
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


%% Parameters cheching

if nargin < 2, errordlg('Incorrect number of input paramters. Please, check the help for further details.'); end
if ~iscell(X), error('The data structure must be a cell array to keep unsynchronized batch trajectories.'); end
    
I =length(X);
l = zeros(1,I);

for i=1:I
    l(i)=size(X{i},1);
end

switch method
    case 'median'
    [num Bref] = min(abs(l - median(l)));
    case 'average'
    [num Bref] = min(abs(l - mean(l)));
    otherwise
        errordlg('The method introduced is not correct.')
end
        