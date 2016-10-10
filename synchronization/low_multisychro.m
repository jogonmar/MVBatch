function [Xs,F,specSynchronization] = low_multisychro(cal,ref,asynDetection,Wconstr,pcs,offset,console,specSynchronization)

% MultiSynchro is devoted to synchronize the key process events ensuring
% the same evolution across batches, no matter the type of asynchronism 
% present in batch data.  
% The multi-synchro algorithm is composed of a high-level and low-level 
% routine. The highlevel % routine is aimed at recognizing the different 
% types of asynchronous trajectories for the subsequent batch classification
% as function of the nature of asynchronism. The low-level routine is in 
% charge of synchronizing the variable trajectories of each one of the 
% batches with a specific procedure based on the type of asynchronism. 
% The original paper are:
%
% [1] González Martínez, JM.; De Noord, O.; Ferrer, A. (2014). 
% Multi-synchro: a novel approach for batch synchronization in scenarios 
% of multiple asynchronisms, Journal of Chemometrics, 28(5):462-475.
%
% [2] González Martínez, JM.; Vitale, R.; De Noord, OE.; Ferrer, A. (2014). 
% Effect of synchronization on bilinear batch process modeling, 
% Industrial and Engineering Chemistry Research, 53(11):4339-4351.
%
% CALLS:
%        [Xs,F,specSynchronization] = MultiSynchro(cal,ref,Wconstr)                                     % minimum call
%        [Xs,F,specSynchronization] = MultiSynchro(cal,ref,Wconstr,offset,console,specSynchronization)  % complete call
%
% INPUTS:
%
% cal: (1xI) cell array containing the measurements collected for J variables at 
%       Ki different sampling times for each one of the I batches.
%       
% ref: (KxJ) reference batch.
%
% Wconstr: (Jx1) boolean array indicating if a specific variables is
% considered in the synchronization (0) or not (1).
%
% pcs: (1x1) number of principal components.
%
% offset: (Jx1) offset to apply to the warping profiles in case that the
% algorithm is applied to different stages (by default array of zeros).
%
% asynDetection: struct containing information derived from the high level routine of the algorithm:
%       - batchcI_II.I: (I1x1) indeces of the batches with class I and II
%       asynchronism.
%       - batchcIII.I: (I2x1) indeces of the batches with class III
%       asynchronism.
%       - batchcIV.I: (I3x1) indeces of the batches with class IV
%       asynchronism.
%       - batchcIII_IV.I: (I3x1) indeces of the batches with class III and IV
%       asynchronism.
% 
% specSynchronization: struct containing the information required to proceed with the
% synchronization of test batches (case that we are monitoring incoming
% batches)
%       - W: (Jx1) non-negative weight array used to give more importance to
%       certain process variables.
%       - Wconstr: (Jx1) boolean array indicating if a specific variable is
%       considered in the synchronization (0) or not (1).
%       - rng: (1xJ) vector containing the mean range of each one the J
%      trajectories.
%       - Xi: (K x J) sample average according to trajectory centering.
%       - Omega: (K x J) sample average according to the trajectory scaling
%       - p: (M x pc) matrix of loadings where pc is the rank of the original X matrix used to fit the PCA model.
%       - t: (N x pc) matrix of scores where pc is the rank of the original X matrix used to fit the PCA model.
%       - offsetnext: (Jx1) offset to apply to the warping profiles in case that the
%       algorithm is applied to another stage.
%
% console: (1x1) handle of the object, 0 for main console.
%
%
% OUTPUTS:
%
% Xs:  (KxJxI) data matrix containing the J batch trajectories, whose duration is equal to
%      K, of each one the I batches.
%
% F: (KrefxI) warping information obtained in the phase of specific synchronization, which is expressed as a function of the
% reference batch.
%
% specSynchronization: struct containing the information required to proceed with the
% synchronization of test batches 
%       - W: (Jx1) non-negative weight array used to give more importance to
%       certain process variables.
%       - Wconstr: (Jx1) boolean array indicating if a specific variable is
%       considered in the synchronization (0) or not (1).
%       - rng: (1xJ) vector containing the mean range of each one the J
%      trajectories.
%       - Xi: (K x J) sample average according to trajectory centering.
%       - Omega: (K x J) sample average according to the trajectory scaling
%       - p: (M x pc) matrix of loadings where pc is the rank of the original X matrix used to fit the PCA model.
%       - t: (N x pc) matrix of scores where pc is the rank of the original X matrix used to fit the PCA model.
%       - offsetnext: (Jx1) offset to apply to the warping profiles in case that the
%       algorithm is applied to another stage.
%
%
% coded by: Jose Maria Gonzalez-Martinez (jogonmar@gmail.com)
%           
% last modification: 21/Aug/14 -> offline and online version of the algorithm are merged. 
%
% Copyright (C) 2014  Technical University of Valencia, Valencia
% Copyright (C) 2014  Jose Maria Gonzalez-Martinez
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
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

if ~iscell(cal), error('The first argument must be a cell array containing the unsynchronized trajectories.'); end
nBatches = size(cal,2);
nTime = size(ref,1);
nVariables = size(ref,2);
if size(ref,2) ~= nVariables, error('The reference batch does not have the same number of variables as the calibration data set.'); end
if ~isstruct(asynDetection), error('asynDetection must be a struct contaning the indices of batches affected by different types of asynchronisms.');end
if ~isfield(asynDetection.batchcI_II,'I'), error('The indices of the batches affected by class I and/or class II asynchronisms.'); end
if ~isfield(asynDetection.batchcIII,'I'), error('The indices of the batches affected by class III asynchronism.'); end
if ~isfield(asynDetection.batchcIV,'I'), error('The indices of the batches affected by class IV asynchronism.'); end
if ~isfield(asynDetection.batchcIII_IV,'I'), error('The indices of the batches affected by class III and class IV asynchronisms.'); end
if nargin < 4, Wconstr = zeros(size(cal{1,1},2),1); end
if nargin < 5 || isempty(pcs), pcs = 6; end
if nargin < 6 || isempty(offset), offset = zeros(nBatches,1);end
if size(offset,1) ~= nBatches, error('Unexpected number offsets. It does not coincide with the number of batches.'); end 
if nargin < 7 || isempty(console), console = 0;end
mode = 1;
if nargin < 8, mode = 0; end
if mode
    if ~isstruct(specSynchronization), error('specSynchronization must be a struct contaning the fields required to proceed with the synchronization of test batches. For further information, type help MultiSynchro on the Matlab prompt');end
    if ~isfield(specSynchronization,'W'), error('The non-negative array W does not exist in the structure specSynchronization.'); end
    if numel(specSynchronization.W) ~= nVariables, error('The number of variables of the batches contained in the test data set does not match with the number of variable weights from the synchronization of the calibration data set.');end
    if ~isfield(specSynchronization,'Wconstr'), error('The array contaning the variable constraints does not exist in the structure specSynchronization.'); end
    if numel(specSynchronization.Wconstr) ~= nVariables, error('The number of variables of the batches contained in the test data set does not match with the number of contraints used in the synchronization of the calibration data set.');end
    if ~isfield(specSynchronization,'rng'), error('The array containg the range values for the process variables rng does not exist in the structure specSynchronization.'); end
    if numel(specSynchronization.rng) ~= nVariables, error('The number of variables of the batches contained in the test data set does not match with the number of ranges calculated in the synchronization of the calibration data set.');end
    if ~isfield(specSynchronization,'Xi'), error('The matrix of averages Xi does not exist in the structure specSynchronization.'); end
    if size(specSynchronization.Xi,2) ~= nVariables, error('The number of variables of the batches contained in the test data set does not match with the number of averages.'); end
    if size(specSynchronization.Xi,1) ~= size(ref,1), error('The number of time points of the reference batch does not coincide with the number of averages in time.'); end
    if ~isfield(specSynchronization,'Omega'), error('The matrix of standard deviations Omega does not exist in the structure specSynchronization.'); end
    if size(specSynchronization.Omega,2) ~= nVariables, error('The number of variables of the batches contained in the test data set does not match with the number of standard deviations.'); end
    if size(specSynchronization.Omega,1) ~= size(ref,1), error('The number of time points of the reference batch does not coincide with the number of standard deviations in time.'); end
    if ~isfield(specSynchronization,'t'), error('The score matrix T does not exist in the structure specSynchronization.'); end
    if ~isfield(specSynchronization,'p'), error('The loading matrix P does not exist in the structure specSynchronization.'); end
else
    if pcs > min(nVariables*nTime,nBatches), error('The number of PCs selected for monitoring exceeds the rank of the batch-wise unfolded matrix');end
end

%% Step C: Specific synchronization

if isempty(asynDetection.batchcI_II.I) && ~mode, errordlg('The Multisynchro approach cannot be applied because no batch with class I and II asynchronism has been found.'); return; end 

    % CLASS I AND II ASYNCHRONISM: DTW-BASED SYNCHRONIZATION/ABNORMALITY DETECTION PROCEDURE

    if ~mode
        cprintMV(console,'Executing the iterative DTW-based synchronization/abnormality detection procedure');
        XG = cal(asynDetection.batchcI_II.I);
        auxbatchcI_II = asynDetection.batchcI_II.I;
        Ift = NaN;

        while(~isempty(Ift))     
            % Step (i): DTW-based synchronization
            [W,X1G,warpX1G,rng] = DTW_Kassidas(XG,ref,Wconstr);
            sGsyn = size(X1G);

            X1G = missTSR3D(X1G,4,sGsyn(1)-1,2);
            % Step (ii): Preprocessing
            [xce, Xi, Omega] = preprocess3D(X1G,2); % Trajectory C&S

            % Step (iii): Batch-wise unfolding and PCA
            xu=unfold(xce,sGsyn(1)-1); % BW unfolding
            pcs = min(pcs,rank(xu));
            [U,S,V] = svd(xu,'econ');
            tAll = U*S;
            pAll = V;
            p = pAll(:,1:pcs);
            t = tAll(:,1:pcs);
            resa=xu-t(:,1:pcs)*p(:,1:pcs)';

            % Step (iv): Design monitoring scheme based on SPE 
            alpha = 0.05;
            m = mean(sum(resa.^2,2));
            v = var(sum(resa.^2,2));
            limr = (v/(2*m))*chi2inv(1-alpha,(2*m^2)/v);

            % Step (v): Off-line post batch monitoring

            % v.1 Compute the SPE statistic for each batch and sort out the
            % corresponding values in ascending order.
            ssqres = sort(sum(resa.^2,2)');

            % v.2 Calculate the acceptable number of batches R that can exceed the control limits at (1-alpha) confidence level by chance as alpha times the number of calibration batches.  
            R = min(1,floor((1-alpha)*sGsyn(3)));

            % v.3 Detect the number of faulty batches in the data set.
            If1 = find(ssqres > limr); if numel(If1 > R), If1 = If1(R+1:end);end
            If2 = find(ssqres > 1.25*limr);
            Ift = union(If1,If2);

            % v.4 Arrange NOC data into the NOC 3-way array
            auxbatchcI_II = setdiff(auxbatchcI_II,auxbatchcI_II(Ift));
            XG = cal(auxbatchcI_II);
            cprintMV(console,'Next iteration of the iterative synchronization');
        end

        % SYNCHRONIZATION OF FAULTY BATCHES FOUND IN THE ITERATIVE MODELING
        indxFaultyBatches = setdiff(asynDetection.batchcI_II.I,auxbatchcI_II);

        XB = cal(indxFaultyBatches);
        X1B = zeros(size(ref,1),size(ref,2),numel(indxFaultyBatches));
        warpXB  = cell(numel(indxFaultyBatches),1);
        warpingXB = cell(numel(indxFaultyBatches),1);
        if ~isempty(indxFaultyBatches)
            for i=1:length(XB)
                [synX1, warpXB{i}, warpingXB{i}]= DTW(scale_(XB{i},rng),scale_(ref,rng),diag(W),false,Wconstr);
                X1B(:,:,i)= synX1.*repmat(rng,size(ref,1),1);
            end
        end

         if numel(asynDetection.batchcIII.I)>0 || numel(asynDetection.batchcIII_IV)>0
             pcs = min(9,min(sGsyn(1),sGsyn(3)));
%         cprintMV(console,'Incompleted batches have been detected. Cross validation is running to determine the optimum number of PCs for missing data imputation.');
%             for a=1:min(50,rank(xu))
%                 cumpress(a) = crossval3D_s(X1G,a,sGsyn(1)-1,ones(sGsyn(1),1));
%             end
%         [~,pcs] = min(cumpress);
         p = pAll(:,1:pcs);
         t = tAll(:,1:pcs);
         end
         
        % STORE INFORMATION OF THE SPECIFIC SYNCHRONIZATION
        specSynchronization.W = W;
        specSynchronization.Wconstr = Wconstr;
        specSynchronization.rng = rng;
        specSynchronization.t = t;
        specSynchronization.p = p;
        specSynchronization.Xi = Xi;
        specSynchronization.Omega = Omega;
        specSynchronization.pcs = pcs;
        
    else
        % CLASS I AND II ASYNCHRONISM: DTW-BASED SYNCHRONIZATION/ABNORMALITY DETECTION PROCEDURE
        if ~isempty(asynDetection.batchcI_II.I),cprintMV(console,'Synchronizing batches with class I and/or II asynchronism.'); end
        X1GB = cal(asynDetection.batchcI_II.I);
        XGB = zeros(size(ref,1),size(ref,2),size(cal,1));
        warpXGB  = cell(numel(asynDetection.batchcI_II.I),1);
        warpingXGB = cell(numel(asynDetection.batchcI_II.I),1);
        if ~isempty(asynDetection.batchcI_II.I)
            for i=1:length(X1GB)
                [synX1, warpXGB{i}, warpingXGB{i}]= DTW(scale_(X1GB{i},specSynchronization.rng),scale_(ref,specSynchronization.rng),diag(specSynchronization.W),false,specSynchronization.Wconstr);
                XGB(:,:,i)= synX1.*repmat(specSynchronization.rng,size(ref,1),1);
            end
        end
    end
    
    % CLASS III AND IV: DTW-BASED SYNCHRONIZATION WITH RELAXED END POINT
    if ~isempty(asynDetection.batchcIII.I),cprintMV(console,'Synchronizing batches with class III asynchronism.'); end
    X2 = cell(numel(asynDetection.batchcIII.I),1);
    warpX2  = cell(numel(asynDetection.batchcIII.I),1);
    warpingX2 = cell(numel(asynDetection.batchcIII.I),1);
    for i=1:numel(asynDetection.batchcIII.I)
        [synX2, warpX2{i}, warpingX2{i}] = DTW(scale_(cal{asynDetection.batchcIII.I(i)},specSynchronization.rng),scale_(ref,specSynchronization.rng),diag(specSynchronization.W),true,specSynchronization.Wconstr);
        X2{i}= synX2.*repmat(specSynchronization.rng,size(synX2,1),1);
    end

    % TSR-based imputation       
     for i=1:size(X2,1)  
        %if  size(X2{i,1},1) < 0.4*nTime
        %     X2s{i} = [X2{i,1}; ones(nTime-size(X2{i,1},1),size(X2{i,1},2)).*NaN];
        %else
        X2s{i} = reconstructX(X2{i,1},specSynchronization.t,specSynchronization.p,specSynchronization.pcs,specSynchronization.Xi,specSynchronization.Omega); 
        %end
     end

    % CLASS IV: DTW-BASED SYNCHRONIZATION WITH RELAXED STARTING POINT
    if ~isempty(asynDetection.batchcIV.I),cprintMV(console,'Synchronizing batches with class IV asynchronism.'); end
    X3 = cell(numel(asynDetection.batchcIV.I),1);
    warpX3  = cell(numel(asynDetection.batchcIV.I),1);
    warpingX3 = cell(numel(asynDetection.batchcIV.I),1);
    for i=1:numel(asynDetection.batchcIV.I)
        [synX3, warpX3{i}, warpingX3{i}] = DTW(scale_(cal{asynDetection.batchcIV.I(i)},specSynchronization.rng),scale_(ref,specSynchronization.rng),diag(specSynchronization.W),[],specSynchronization.Wconstr);
        X3{i}= synX3.*repmat(specSynchronization.rng,size(synX3,1),1);
    end

    % CLASS III AND IV: DTW-BASED SYNCHRONIZATION WITH RELAXED STARTING AND END POINT
    if ~isempty(asynDetection.batchcIII_IV.I),cprintMV(console,'Synchronizing batches with class III and IV asynchronisms.'); end
    X4 = cell(numel(asynDetection.batchcIII_IV.I),1);
    warpX4  = cell(numel(asynDetection.batchcIII_IV.I),1);
    warpingX4 = cell(numel(asynDetection.batchcIII_IV.I),1);
    for i=1:numel(asynDetection.batchcIII_IV.I)
        [synX4, warpX4{i}, warpingX4{i}] = DTW(scale_(cal{asynDetection.batchcIII_IV.I(i)},specSynchronization.rng),scale_(ref,specSynchronization.rng),diag(specSynchronization.W),true,specSynchronization.Wconstr);
        X4{i}= synX4.*repmat(specSynchronization.rng,size(synX4,1),1);
    end

    % TSR-based imputation of CLASS III AND IV aynchronisms
%      for i=1:size(X4,1)   
%         X4s{i} = reconstructX(X4{i,1},specSynchronization.t,specSynchronization.p,specSynchronization.pcs,specSynchronization.Xi,specSynchronization.Omega); 
%      end
     
      for i=1:size(X4,1)  
        %if  size(X4{i,1},1) < 0.4*nTime
        %     X4s{i} = [X4{i,1}; ones(nTime-size(X4{i,1},1),size(X4{i,1},2)).*NaN];
        %else
        X4s{i} = reconstructX(X4{i,1},specSynchronization.t,specSynchronization.p,specSynchronization.pcs,specSynchronization.Xi,specSynchronization.Omega); 
        %end
     end

%% MERGE STEP 
 
Xs = zeros(size(specSynchronization.Xi,1),size(specSynchronization.Xi,2),nBatches);
F = ones(size(specSynchronization.Xi,1),nBatches).*NaN;

offsetnext = zeros(nBatches,1);

for i=1:nBatches
    if ~mode
        if ~isempty(find(auxbatchcI_II==i)), Xs(:,:,i) = X1G(:,:,find(auxbatchcI_II==i)); F(:,i) = warpX1G(:,find(auxbatchcI_II==i))+ones(size(warpX1G,1),1).*offset(i);offsetnext(i)=F(end,i); end 
        if ~isempty(find(indxFaultyBatches==i)), Xs(:,:,i) = X1B(:,:,find(indxFaultyBatches==i)); F(:,i) = warpXB{find(indxFaultyBatches==i)}+ones(size(warpXB{find(indxFaultyBatches==i)},1),1).*offset(i); offsetnext(i)=F(end,i); end
    else
        if ~isempty(find(asynDetection.batchcI_II.I==i)), Xs(:,:,i) = XGB(:,:,find(asynDetection.batchcI_II.I==i)); F(:,i) = warpXGB{find(asynDetection.batchcI_II.I==i)}+ones(size(warpXGB{find(asynDetection.batchcI_II.I==i)},1),1).*offset(i); offsetnext(i)=F(end,i); end   
    end
    if ~isempty(find(asynDetection.batchcIII.I==i)), Xs(:,:,i) = X2s{find(asynDetection.batchcIII.I==i)};F(:,i) = warpX2{find(asynDetection.batchcIII.I==i)}+ones(size(warpX2{find(asynDetection.batchcIII.I==i)},1),1).*offset(i); offsetnext(i)=F(end,i); end
    if ~isempty(find(asynDetection.batchcIV.I==i)), Xs(:,:,i) = X3{find(asynDetection.batchcIV.I==i)};F(:,i) = warpX3{find(asynDetection.batchcIV.I==i)}+ones(size(warpX3{find(asynDetection.batchcIV.I==i)},1),1).*offset(i); offsetnext(i)=F(end,i); end
    if ~isempty(find(asynDetection.batchcIII_IV.I==i)), Xs(:,:,i) = X4s{find(asynDetection.batchcIII_IV.I==i)};F(:,i) = warpX4{find(asynDetection.batchcIII_IV.I==i)}+ones(size(warpX4{find(asynDetection.batchcIII_IV.I==i)},1),1).*offset(i); offsetnext(i)=F(end,i); end
end

% Store offset for possible next stage synchronization
specSynchronization.offset = offsetnext;
% Estimate the boundaries of the warping information
specSynchronization.band = estimationBD(F);