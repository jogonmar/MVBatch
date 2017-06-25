classdef LatentModel
% Class for the latent structure used for the design of a monitoring system
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
        % Model parameters
        T          % (phases x 1) cell array containing the score matrices of each model, one per phase.
        P          % (phases x 1) cell array containing the loading matrices of each model, one per phase.
        prep       % prep: (1x1) preprocesing of the data
                   %       0: no preprocessing.
                   %       1: trajectory centering (average trajectory subtraction)
                   %       2: 1 + trajectory-scaling (scales data so that each pair variable and 
                   %           sampling time has variance 1) (default)  
                   %       3: 1 + variable-scaling (scales data so that each variable has
                   %           variance 1)
                   %       4: variable centering (subtraction of the average value of each
                   %           variable)
                   %       5: 4 + variable-scaling.   
        mn         % (K x J) sample average according to the preprocessing method.
        stnd       % (K x J) sample scale according to the preprocessing method.
        phases     % (1x5) phases of the MP model: [PRESS, PCs, lags, initial time, end time].
        dimensions % structure with three fields:
                     % nbatches   (1x1): number of batches of the calibration data set.
                     % nvariables (1x1): number of variables of the calibration data set.
                     % ntimes     (1x1): number of time points of the calibration data set.
        % Monitoring parameters
        monitoring_statistics  % instance of the class Monitoring Parameters
        control_limits         % instance of the class ControlLimits
        % Cross-validated monitoring parameters
        cvD          % (Ix1) cross-validated overall D-statistic values
        cvQ          % (Ix1) cross-validated overall Q-statistic values 
        cvevolD      % (KxI) cross-validated online D-statistic values for I batches
        cvevolQ      % (KxI) cross-validated online Q-statistic values for I batches
    end
    
    methods (Access = public)      
        
        function class_object = LatentModel(xini,phases)
        % Constructor function of the class 
        %
        % INPUTS:
        %
        % xini: (KxJxI) three-way batch data matrix for calibration, K(sampling times) 
        %       x J(variables) x I(batches)
        %
        % phases: (n_phasesx5) phases of the MP model. Each row contains the information 
        %   of a phase, namely [PRESS, PCs, lags, initial time, end time]. 
        %
        %
        % OUTPUTS:
        %
        % class_object: (1x1) instance of the current class.

        % Parameter checking 
        assert(nargin == 2, 'Incorrect number of input parameters');
        assert (ndims(xini)==3, 'Incorrect number of dimensions of xini.');
        s = size(xini);
        assert (isempty(find(s<1)), 'Incorrect number of dimensions of xini.');
        assert(length(s) == 3, 'Incorrect number of parameters of the 3-way array dimension');
        assert(ndims(phases)==2,'Incorrect number of dimensions of phases.');
        sp=size(phases);
        assert(sp(1)>0 && sp(2)==5,'Incorrect content of phases.');
        assert(isempty(find(phases(:,1:3)<0)),'Incorrect content of phases.');
        assert(isempty(find(phases(:,4:5)<1)),'Incorrect content of phases.');
        assert(isempty(find(phases(:,3:5)>s(1))),'Incorrect content of phases.');

        % Model meta-parameters
        class_object.T = cell(sp(1),1);       
        class_object.P = cell(sp(1),1);  
        class_object.prep = 2;
        class_object.mn = cell(sp(1),1);         
        class_object.stnd = cell(sp(1),1);     
        class_object.dimensions = struct('nbatches',s(3),'ntimes',s(1),'nvariables',s(2));
        class_object.phases = phases;
        % Monitoring parameters
        class_object.monitoring_statistics = LatentStructure.MonitoringParameters('calibration');
        class_object.control_limits = LatentStructure.ControlLimits;
        % Cross-validated monitoring parameters
        class_object.cvD = [];        
        class_object.cvQ = [];          
        class_object.cvevolD = [];      
        class_object.cvevolQ = [];
        end
    end
end