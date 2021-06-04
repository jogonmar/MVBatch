function [sys,x0,str,ts] = sfunc(t,x,u,flag,xini, k1h, K1h, k1l, K1l, k1e, K1e, K1i, k2, K2, K2i, k3, K3, k4, K4, k5, K5, k5e, K5e, K5i, k6, K6, k6r, K6e, k7, K7, k8, k9, K9, k9e, K9e, K9i, k9c, k10, K10, k10e, K10e, k11, Sf, w, var_ini)

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


if nargin==0, 
    xini([1 6:8]) = var_ini/(2*100)*randn(1).*xini([1 6:8]) + xini([1 6:8]); % variabilidad inicial
    sys=[8,0,10,1,0,0,1];
    x0=xini;
    str=[];
    ts=[0 0];
    return, 
end

if abs(flag) == 0,     
    xini([1 6:8]) = var_ini/(2*100)*randn(1).*xini([1 6:8]) + xini([1 6:8]); % variabilidad inicial
    sys=[8,0,10,1,0,0,1];
    x0=xini;
    str=[];
    ts=[0 0];
    
else
    
    
    for i=1:8,
        if x(i) < 0, x(i) = 0; end;
    end;
    
    % Model reactions
    ind_t = min(max(ceil(t),1),35);
    
    r1 = k1l(ind_t)*(x(1)/(x(1)+K1l(ind_t)))*x(7) + k1h(ind_t)*(x(1)/(x(1)+K1h(ind_t)))*x(7) + k1e(ind_t)*(x(1)/(x(1)*(K1i(ind_t)*x(3)+1)+K1e(ind_t)))*x(3)*x(7);
    r2 = k2(ind_t)*(x(2)/(x(2)+K2(ind_t)))*(1/(K2i(ind_t)*x(1)+1))*x(7);
    r3 = k3(ind_t)*((x(2)^4)/((x(2)^4)+K3(ind_t)))*x(7);
    r4 = k4(ind_t)*(x(3)/(x(3)+K4(ind_t)))*x(7)*x(8);
    r5 = k5(ind_t)*(x(4)/(x(4)+K5(ind_t)))*x(7) + k5e(ind_t)*(x(4)/(x(4)+K5e(ind_t)))*(1/(1+K5i(ind_t)*x(1)))*x(7);
    r6 = k6(ind_t)*((x(3)-k6r(ind_t)*x(5))/(x(3)+K6(ind_t)+K6e(ind_t)*x(5)))*x(7);
    r7 = k7(ind_t)*(x(1)/(x(1)+K7(ind_t)))*x(7);
    r8 = k8(ind_t)*(x(4)/(x(4)+K5e(ind_t)))*(1/(1+K5i(ind_t)*x(1))*x(7));
    r9 = (k9(ind_t)*(x(1)/(x(1)+K9(ind_t))) + k9e(ind_t)*(x(5)/(x(5)+K9e(ind_t))) )*(1/(K9i(ind_t)*x(1)+1))*x(7) + k9c(ind_t)*(x(1)/(x(1)+K9(ind_t)))*x(7);
    
    r10 = k10(ind_t)*(x(1)/(x(1)+K10(ind_t)))*x(7) + k10e(ind_t)*(x(5)/(x(5)+K10e(ind_t)))*x(7);
    r11 = k11(ind_t)*x(8);
    mu = 0.732*r7 + 0.619*r8;
    
    D = 0;
    
if abs(flag) == 1,
        
    % Mass balances
    sys(1) = -(r1 + r7)*x(6) + (Sf-x(1))*D;
    sys(2) = (0.978*r1 - r2 - r3)*x(6) - x(2)*D;
    sys(3) = (0.5*r3 - r4 - r6)*x(6) - x(3)*D;
    sys(4) = (1.363*r4 - r5 - r8)*x(6) - x(4)*D;
    sys(5) = 1.045*r6*x(6) - x(5)*D;
    sys(6) = (mu - D)*x(6);        
    sys(7) = 0.732*r7 + 0.619*r8 - r9 - r10 - (0.732*r7 + 0.619*r8)*x(7);
    sys(8) = r9 - r11 - (0.732*r7 + 0.619*r8)*x(8);
    
    
elseif flag == 3,
    
    % Salidas
    sys(1:8) = x(1:8);
    sys(9) = (1000/32)*(.178*r1 + .908*r2 + .363*r4 + 1.066*r5 - 0.363*r6 + 0.063*r7 + .214*r8);
    sys(10) = (1000/44.01)*(1.499*r2 + .5*r3 + 1.466*r5 + .127*r7 + .325*r8);

else 
    sys = []; 
end
end