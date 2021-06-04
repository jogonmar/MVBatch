% Batch process simulation: Fermentation of Saccharomyces cerevisiae
%
% Model used: Lei F., M. Rotboll, Jorgensen S.B. A biochemically structured
%   model for Saccharomyces cerevisiae. Journal of Biotechnology. 2001.
%   88:205-221.
%
% Process Variables:
%
%   1 - Glucose concentration
%   2 - Pyruvate concentration
%   3 - Acetald concentration
%   4 - Acetate concentration
%   5 - Ethanol concentration
%   6 - Biomass concentration (dry weight)
%   7 - Active cell material
%   8 - Acetaldehyde dehydrogenase proportional to the measured activity
%   9 - Specific oxygen uptake rate
%   10 - Specific carbon dioxide evolution rate
%
%
% Citation:
% J.M. González-Martínez, J. Camacho, and A. Ferrer. MVbatch: a matlab toolbox for batch process modeling and monitoring. Chemometrics and Intelligent Laboratory Systems, 183:122-133, 2018.
%
% Copyright (C) 2021  Jose Camacho Paez, University of Granada, Granada
% Copyright (C) 2021  José M. González Martínez
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


clear
n_cal=30;   % number of calibration batches
n_tNOC=1;  % number of test batches under NOC
n_tab=1;    % number of test abnormal batches, there are 6 different abnormalities
ini         % parameters of the model
% seed to control the simulation, change this value for a different
% simulation
stream = RandStream('mt19937ar','seed',9);
RandStream.setGlobalStream(stream); 

% Initial state
xini = [30,0,0,0,0,.002,0.1,.0075];

var_ini = 10;   % 10% of initial variability
med_err = 5;    % 5% of measurement error (all outputs) 


%% Calibration and NOC test bacthes

for i=1:n_cal,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    cal{i}=[noise_corrupted.signals.values noise_corrupted.time]
end

for i=1:n_tNOC,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    testg{i}=[noise_corrupted.signals.values noise_corrupted.time]
end


%% Abnormalities 1 and 2: change in the k1l constant

k1l=2.43*ones(35,1);

for i=1:n_tab,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    testb1{i}=[noise_corrupted.signals.values noise_corrupted.time]
end

k1l=0.43*ones(35,1);

for i=1:n_tab,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    testb2{i}=[noise_corrupted.signals.values noise_corrupted.time]
end

k1l=1.43*ones(35,1);


%% Abnormalities 3 and 4: change in the k6 constant

k6=[3.8*ones(15,1);2.2*ones(20,1)];

for i=1:n_tab,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    testb3{i}=[noise_corrupted.signals.values noise_corrupted.time]
end


k6=[2.32*ones(15,1);2.82*ones(20,1)];

for i=1:n_tab,
    r=round(abs(100*randn(10,1)));
    sim('saccha');
    testb4{i}=[noise_corrupted.signals.values noise_corrupted.time]
end

k6=2.82*ones(35,1);


%% Abnormalities 5 and 6: sensor bias in the biomass concentration

testb5 = testg(1:n_tab);
testb6 = testg(1:n_tab);

for i=1:n_tab,
    indx = find(testb5{i}(:,end)>10 & testb5{i}(:,end)<=20);
    testb5{i}(indx,6) = testb5{i}(indx,6)+1;
    testb6{i}(indx,6) = testb6{i}(indx,6)-1;
end


%% Obtain data with non-equalized variables

unequalize

%% Save simulated data

save saccha