classdef ControlLimits 
% Class for the monitoring parameters associated with a monitoring system
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
        % suggested and cv imposed significance level (alpha) for post-batch offline monitoring
        % alpd95: (1x1) 95% confidence level for the overall D-statistic. 
        % alpq95: (1x1) 95% confidence level for the overall Q-statistic.
        % alpd99: (1x1) 99% confidence level for the overall D-statistic. 
        % alpq99: (1x1) 99% confidence level for the overall Q-statistic.
        % alpd95cv: (1x1) 95% confidence level for the cross-validated overall D-statistic. 
        % alpq95cv: (1x1) 95% confidence level for the cross-validated overall Q-statistic.
        % alpd99cv: (1x1) 99% confidence level for the cross-validated overall D-statistic. 
        % alpq99cv: (1x1) 99% confidence level for the cross-validated overall Q-statistic.
        offline_alpha = struct('alpd95',0.05,'alpq95',0.05,'alpd99',0.01,'alpq99',0.01,...
                          'alpd95cv',Inf,'alpq95cv',Inf,'alpd99cv',Inf,'alpq99cv',Inf);
                      
        % suggested and cv imposed significance level (alpha) for post-batch online monitoring
        % alpod95: (1x1) 95% confidence level for the online D-statistic. 
        % alpoq95: (1x1) 95% confidence level for the online Q-statistic. 
        % alpod99: (1x1) 99% confidence level for the online D-statistic. 
        % alpoq99: (1x1) 99% confidence level for the online Q-statistic.
        % alpod95cv: (1x1) 95% confidence level for the cross-validated online D-statistic. 
        % alpoq95cv: (1x1) 95% confidence level for the cross-validated online Q-statistic. 
        % alpod99cv: (1x1) 99% confidence level for the cross-validated online D-statistic. 
        % alpoq99cv: (1x1) 99% confidence level for the cross-validated online Q-statistic.
        online_alpha  = struct('alpd95',0.05,'alpq95',0.05,'alpd99',0.01,'alpq99',0.01,...
                          'alpd95cv',Inf,'alpq95cv',Inf,'alpd99cv',Inf,'alpq99cv',Inf);     
        
        % theoretical and cv control limits for post-batch offline monitoring 
        % limd95: (1x1) 95% confidence limit for the overall D-statistic. 
        % limq95: (1x1) 95% confidence limit for the overall Q-statistic.
        % limd99: (1x1) 99% confidence limit for the overall D-statistic. 
        % limq99: (1x1) 99% confidence limit for the overall Q-statistic.
        % limd95cv: (1x1) 95% confidence limit for the cross-validated overall D-statistic. 
        % limq95cv: (1x1) 95% confidence limit for the cross-validated overall Q-statistic.
        % limd99cv: (1x1) 99% confidence limit for the cross-validated overall D-statistic. 
        % limq99cv: (1x1) 99% confidence limit for the cross-validated overall Q-statistic.
        offline_cl = struct('limd95',[],'limq95',[],'limd99',[],'limq99',[],...
                  'limd95cv',[],'limq95cv',[],'limd99cv',[],'limq99cv',[]);              
     
        % theoretical and cv control limits for post-batch online monitoring
        % limod95: (1x1) 95% confidence limit for the online D-statistic. 
        % limoq95: (1x1) 95% confidence limit for the online Q-statistic.
        % limod99: (1x1) 99% confidence limit for the online D-statistic. 
        % limoq99: (1x1) 99% confidence limit for the online Q-statistic. 
        % limod95cv: (1x1) 95% confidence limit for the cross-validated online D-statistic. 
        % limoq95cv: (1x1) 95% confidence limit for the cross-validated online Q-statistic.
        % limod99cv: (1x1) 99% confidence limit for the cross-validated online D-statistic. 
        % limoq99cv: (1x1) 99% confidence limit for the cross-validated online Q-statistic.
        online_cl = struct('limd95',[],'limq95',[],'limd99',[],'limq99',[],...
                  'limd95cv',[],'limq95cv',[],'limd99cv',[],'limq99cv',[]);       
    end
end