% Batch process simulation: Fermentation of Saccharomyces cerevisiae
%
% Model used: Lei F., M. Rotboll, Jorgensen S.B. A biochemically structured
%   model for Saccharomyces cerevisiae. Journal of Biotechnology. 2001.
%   88:205-221.
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


% Constants
k1h = .584*ones(35,1);
K1h = .0116*ones(35,1);
k1l = 1.43*ones(35,1);
K1l = .94*ones(35,1);
k1e = 47.1*ones(35,1);
K1e = .12*ones(35,1);
K1i = 14.2*ones(35,1);
k2 = .501*ones(35,1);
K2 = 2e-5*ones(35,1);
K2i = .101*ones(35,1);
k3 = 5.81*ones(35,1);
K3 = 5e-7*ones(35,1);
k4 = 4.8*ones(35,1);
K4 = 2.64e-4*ones(35,1);
k5 = .0104*ones(35,1);
K5 = .0102*ones(35,1);
k5e = .775*ones(35,1);
K5e = .1*ones(35,1);
K5i = 440*ones(35,1);
k6 = 2.82*ones(35,1);
K6 = .034*ones(35,1);
k6r = .0125*ones(35,1);
K6e = .057*ones(35,1);
k7 = 1.203*ones(35,1);
K7 = .0101*ones(35,1);
k8 = .589*ones(35,1);
k9 = .008*ones(35,1);
K9 = 1e-6*ones(35,1);
k9e = .0751*ones(35,1);
K9e = 13*ones(35,1);
K9i = 25*ones(35,1);
k9c = 3.99e-3*ones(35,1);
k10 = .392*ones(35,1);
K10 = 2.3e-3*ones(35,1);
k10e = 3.39e-3*ones(35,1);
K10e = 1.8e-3*ones(35,1);
k11 = .02*ones(35,1);

Sf = 75;
w=3;
