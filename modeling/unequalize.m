% Batch process simulation: Fermentation of Saccharomyces cerevisiae
% Obtain data with non-equalized variables
%
% Model used: Lei F., M. Rotboll, Jorgensen S.B. A biochemically structured
%   model for Saccharomyces cerevisiae. Journal of Biotechnology. 2001.
%   88:205-221.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/Aug/09
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

cal_neq = [];
for i=1:length(cal),
    t = 1:size(cal{i},1);
    calstr.data{1} = [t(1:2:end)' cal{i}(1:2:end,1:3)];
    calstr.data{2} = [t(2:3:end)' cal{i}(2:3:end,4:6)];
    calstr.data{3} = [t(3:5:end)' cal{i}(3:5:end,7:11)];
    cal_neq = [cal_neq calstr];
end

testg_neq = [];
for i=1:length(testg),
    t = 1:size(testg{i},1);
    teststr.data{1} = [t(1:2:end)' testg{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testg{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testg{i}(3:5:end,7:11)];
    testg_neq = [testg_neq teststr];
end

testb1_neq = [];
for i=1:length(testb1),
    t = 1:size(testb1{i},1);
    teststr.data{1} = [t(1:2:end)' testb1{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb1{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb1{i}(3:5:end,7:11)];
    testb1_neq = [testb1_neq teststr];
end

testb2_neq = [];
for i=1:length(testb2),
    t = 1:size(testb2{i},1);
    teststr.data{1} = [t(1:2:end)' testb2{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb2{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb2{i}(3:5:end,7:11)];
    testb2_neq = [testb2_neq teststr];
end

testb3_neq = [];
for i=1:length(testb3),
    t = 1:size(testb3{i},1);
    teststr.data{1} = [t(1:2:end)' testb3{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb3{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb3{i}(3:5:end,7:11)];
    testb3_neq = [testb3_neq teststr];
end

testb4_neq = [];
for i=1:length(testb4),
    t = 1:size(testb4{i},1);
    teststr.data{1} = [t(1:2:end)' testb4{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb4{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb4{i}(3:5:end,7:11)];
    testb4_neq = [testb4_neq teststr];
end

testb5_neq = [];
for i=1:length(testb5),
    t = 1:size(testb5{i},1);
    teststr.data{1} = [t(1:2:end)' testb5{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb5{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb5{i}(3:5:end,7:11)];
    testb5_neq = [testb5_neq teststr];
end

testb6_neq = [];
for i=1:length(testb6),
    t = 1:size(testb6{i},1);
    teststr.data{1} = [t(1:2:end)' testb6{i}(1:2:end,1:3)];
    teststr.data{2} = [t(2:3:end)' testb6{i}(2:3:end,4:6)];
    teststr.data{3} = [t(3:5:end)' testb6{i}(3:5:end,7:11)];
    testb6_neq = [testb6_neq teststr];
end
