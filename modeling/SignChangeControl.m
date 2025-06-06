function [p_corr,p2_corr] = SignChangeControl(p,p2)

% Sign correction to account for potential rotational ambiguity inherent in 
% PCA models so that loadings that can be compared across models. 
%
% [p,p2] = SignChangeControl(p,p2) % minimum call
%
% INPUTS:
%
% p:  [MxA] loading matrix corresponding to first latent-variable model
%
% p2: [MxA] loading matrix corresponding to first latent-variable model
%
%
% OUTPUTS:
%
% p_corr:  [MxA] sign-corrected loading matrix corresponding to first latent-variable model
%
% p2_corr: [MxA] sign-corrected loading matrix corresponding to first latent-variable model
%
%
% coded by: Jose M. Gonzalez Martinez (jogonmar@gmail.com)
% last modification: 10/Jan/2025
%
% Copyright (C) 2014  José M. Gonzalez-Martinez
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


%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin == 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[J_p,A_p]   = size(p);
[J_p2,A_p2] = size(p2);
assert (J_p == J_p2, 'Error: Loading matrices belong to latent-variable models fitted on two-way arrays containing different number of components. Type ''help %s'' for more info.', routine(1).name);
assert (A_p == A_p2, 'Error: Loading matrices contain different number of components. Type ''help %s'' for more info.', routine(1).name);

[v ind] =  max(abs(p));
v2= p2(ind);

sg = sign(p(ind));
sg2 = sign(p2(ind));

if sg~=sg2
    if sg==1, sg2=-1;
    else
        sg=1;
        sg2=-1;
    end
else
    sg=1;sg2=1;
end

p_corr = sg.*p;
p2_corr = sg2.*p2;