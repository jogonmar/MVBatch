function [calibration] = BuildDataBase(cal,varNames,batchNames)

if ~iscell(cal), error('The calibration dataset must be provided in a cell structure'); end
nVariables = size(cal{1,1},2);
if nargin < 2
    nVariables = nVariables + 2;
    varNames = cell(nVariables,1);
    varNames{1,1} = 'Time';
    varNames{2,1} = 'Stage';
    for i=3:nVariables
        varNames{i,1} = strcat('Variable ',num2str(i-2));
    end
end
if length(varNames(:,1))~= nVariables, error('The number of tagnames does not match with the number of variables of the dataset'); end

nBatches = length(cal);
if nargin < 3
    
    batchNames = cell(nBatches,1);
    for i=1:nBatches
        batchNames{i,1} = num2str(i);
    end
end
if length(batchNames)~=nBatches, error('The number of batch IDs does not match with the number of batches of the dataset');end
calibration.batch_data = struct('data',{});

for i=1:length(cal)
    calibration.batch_data(i,1).data{1,1} = [[1:size(cal{i},1)]' ones(size(cal{i},1),1) cal{i}];
end
calibration.varNames = [varNames repmat({'Batch time','Units'},nVariables,1)];
calibration.batchNames = batchNames;