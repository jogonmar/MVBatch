% Batch process simulation: Fermentation of Saccharomyces cerevisiae
% Obtain data with non-equalized variables
%
% Model used: Lei F., M. Rotboll, Jorgensen S.B. A biochemically structured
%   model for Saccharomyces cerevisiae. Journal of Biotechnology. 2001.
%   88:205-221.
%
% codified by: José Camacho Páez.
% version: 0.0
% last modification: 21/Aug/09.

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
