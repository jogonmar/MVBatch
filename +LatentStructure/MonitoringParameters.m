classdef MonitoringParameters
% Class for the multivariate statistics of calibration data sets
% coded by: José M. González Martínez (jogonmar@gmail.com)
% Copyright (C) 2017  José M. González Martínez
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
    
    properties 
        % Monitoring parameters
        D           = [];         % (Ix1): overall D-statistic values for I test batches
        Q           = [];         % (Ix1): overall Q-statistic values for I test batches 
        evolD       = [];         % (KxJ): K online D-statistic values for I test batches
        evolQ       = [];         % (KxJ): K online D-statistic values for I test batches 
        cont_evolD  = [];     	  % (KxIJ): (KxJ) contributions to the D statistic for I test batches 
        cont_evolQ  = [];         % (KxIJ): (KxJ) contributions to the Q statistic for I test batches 
        name        = '';         %  Name of the test data set
    end
    
    methods (Access = public)      
        
        function class_object = MonitoringParameters(name)
            % Constructor function of the class 
            % Monitoring parameters
            class_object.D = [];        
            class_object.Q = [];          
            class_object.evolD = [];      
            class_object.evolQ = [];        
            class_object.cont_evolD = [];
            class_object.cont_evolQ = [];
            class_object.cont_evolQ = [];
            class_object.name = name;
        end
        
    end
end