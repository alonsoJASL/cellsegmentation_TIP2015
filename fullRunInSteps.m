function fullRunInSteps(imCytoSet, imNucleiSet, storageCommonPath, ...
                        inparameters, outputpath)
% Run the function in steps.
 
%%  STEP 1 - Cytoplasm's (clump) Mask by Superpixel & Convex Hull
try
    load(strcat(outputpath, storageCommonPath, 'ConvexhullCytoClumpMaskSet.mat'), ...
        'ConvexhullCytoClumpMaskSet');
    fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement.........done!\n');
catch
    fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement...\n');
    
    ConvexhullCytoClumpMaskSet = cell(imNum,1);

    % Superpixel Parameters by Data set
    spRATIO = 0.5;
    spKERNELSIZE = 2;
    spMAXDIST = 10;
    
    % Timer
    t_preproc = zeros(size(imCytoSet,1),1);
    for i = 1:length(imCytoSet)
        tic;
        im_i = imCytoSet{i,1};
        
        fprintf('\tSuperpixel Convex Hull: Image %d out of %d \n',...
        i, imNum);
        [ConvexhullCytoClumpMask] = SuperpixelConvexHull( im_i, ...
            spRATIO, spKERNELSIZE, spMAXDIST );
        ConvexhullCytoClumpMaskSet{i,1} = logical(ConvexhullCytoClumpMask);
        
        t_preproc(i) = toc; 
    end
    
    fprintf('\tStep - 1 ...... DONE! \n');
    %| Save section variable to file |
    save(strcat(outputpath, storageCommonPath, 'ConvexhullCytoClumpMaskSet.mat'),...
        'ConvexhullCytoClumpMaskSet', 't_preproc');
end
%% STEP 2 - Cytoplasm's (clump) Mask by GMM Learning
try
    load(strcat(outputpath, storageCommonPath, 'GMMCytoClumpMaskSet.mat'), ...
        'GMMCytoClumpMaskSet', 'gmm_model_cytoClump', ...
        'gmm_model_cytobackground', 't_GMM_train');
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary.........done!\n');
catch
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary...\n');
   
    % Timer 
    tic;
    
    % Semi-supervised Learning by GMM
    idx_TrainImg = (1:imNum)';
    GMMCytoClumpMaskSet = GenerateForeBackgroundMask(ConvexhullCytoClumpMaskSet);
    
    for ixGMM = 1:10
        if ixGMM > 1
            GMMCytoClumpMaskSet = GenerateForeBackgroundMask( GMMCytoClumpMaskSet );
        end
        
        % Train GMM model by ALL images
        [gmm_model_cytoClump, gmm_model_cytobackground]  = ...
            TrainGMMSceneSegmentation(imCytoSet, GMMCytoClumpMaskSet, ...
                                    idx_TrainImg); 
                                
        % Test GMM model by ALL images
        for i = 1:imNum
            im_i = imCytoSet{i,1};  
            % GMM Posterior for Each Pixel
            [gmm_post] = TestGMMSceneSegmentation(gmm_model_cytoClump,...
                gmm_model_cytobackground, imCytoSet, i);
            
            [GMMCytoClumpMask] = ComputeConfidenceAsScene(gmm_post{1,1}, ...
                size(im_i,1), size(im_i,2));
            GMMCytoClumpMaskSet{i,1} = GMMCytoClumpMask; 
        end        
        fprintf('\tGMM Learning Iteration: %d\n', ixGMM);
    end
   
    t_GMMtrain = toc;  
    
    % Save section variable to file
    save(strcat(outputpath, storageCommonPath, 'GMMCytoClumpMaskSet.mat'), ...
        'GMMCytoClumpMaskSet', 'gmm_model_cytoClump', ...
        'gmm_model_cytobackground', 't_GMMtrain');
        
end

%% STEP 2.1 - Remove Fragments in GMM Learning Result

% Level set paramters for clump refinement (after GMM)
iter_in = inparameters.iter_in;
iter_out = inparameters.iter_out;
alfaGMM = inparameters.alfaGMM;
lambdaGMM = inparameters.lambdaGMM;

try
    load(strcat(outputpath, storageCommonPath, 'SceneCytoClumpMaskSet.mat'), ...
        'SceneCytoClumpMaskSet', 't_GMM_LevelSetClean');
    fprintf('Step - 2.1 Clean Noise in GMM Learned Cyto Clump Masks.........done!\n');
catch
    fprintf('Step - 2.1 Clean Noise in GMM Learned Cyto Clump Masks...\n');
    
    % Timer
    t_gmmclean = zeros(imNum,1);
    
    % Remove Fragments by Level Set Method
    SceneCytoClumpMaskSet = cell(imNum,1);
    for i = 1:imNum   
        tic;
        im_i = imCytoSet{i,1};
        GMMCytoClumpMask = GMMCytoClumpMaskSet{i,1};
        
        % Erase Tiny Fragments by Level Set Method
        GMMCytoClumpMask = cleanFragmentsLevelSet(...
                                            im_i, GMMCytoClumpMask, ...
                                            iter_in, iter_out, alfaGMM, lambdaGMM );
        GMMCytoClumpMask = im2bw(GMMCytoClumpMask);        
        imCBMaskLabels = bwlabel(~GMMCytoClumpMask,8);
        
        % Erase Still Existing Tiny Fragments by Threshold
        for j = 1:max(imCBMaskLabels(:))
            idx = find(imCBMaskLabels == j);
            if size(idx,1) < inparameters.minCellSize
                GMMCytoClumpMask(idx) = 1;
            end
        end
        GMMCytoClumpMask = im2bw(GMMCytoClumpMask);
        
        % Refined Cyto Clump Mask 
        SceneCytoClumpMaskSet{i,1} = ~GMMCytoClumpMask;
 
        t_gmmclean(i) = toc;   
    end
    
    % Save section variable to file
    save(strcat(outputpath, storageCommonPath, 'SceneCytoClumpMaskSet.mat'), ...
        'SceneCytoClumpMaskSet', 't_gmmclean');
    fprintf('\tStep - 2 ...... DONE! \n');
end
 
%% Using MSER to find Nuclei Blobs
%==================================================================

% Level set parameters for nuclei refinement.
iter_in_rawNuclei = inparameters.iter_in_rawNuclei;
iter_out_rawNuclei = inparameters.iter_out_rawNuclei;
alfa_rawNuclei = inparameters.alfa_rawNuclei;
lambda_rawNuclei = inparameters.lambda_rawNuclei;

try
    load(strcat(outputpath, storageCommonPath, 'MSERNucleiBlobSet.mat'),....
        'MSERNucleiBlobSet');
    disp('Step - 3 To find nuclei.........done!');
catch
    disp('Step - 3 To find nuclei...');
    
    % MSER Parameters by Data set
    RegionAreaRange = [120,600];
    ThresholdDelta = 3;

    % MSER Blobs Set 
    MSERNucleiBlobSet = cell(length(imNucleiSet),1);
 
    % timer
    t_mserblobs = zeros(size(imNucleiSet,1),1);
    
    for i = 1:length(imNucleiSet)
        tic;
        I = imNucleiSet{i,1};
 
        % Using MSER to find raw candidate regions for nuclei
        mser_regions = detectMSERFeatures(I, 'RegionAreaRange', ...
            RegionAreaRange, 'ThresholdDelta', ThresholdDelta);
        MSERNucleiBlob = false(size(I));
        for j = 1:length(mser_regions)    
            mser_region_xy = mser_regions(j,1).PixelList;
            for k = 1:size(mser_regions(j,1).PixelList,1)
                IIx = mser_region_xy(k,1);
                IIy = mser_region_xy(k,2);
                MSERNucleiBlob(IIy,IIx) = 1;
            end
        end
 
        MSERNucleiBlob = logical(MSERNucleiBlob);
 
        % Denoise by Level Set Method 
        MSERNucleiBlob = cleanFragmentsLevelSet(I, MSERNucleiBlob, ... 
                                                iter_in_rawNuclei, iter_out_rawNuclei, ...
                                                alfa_rawNuclei, lambda_rawNuclei );

        MSERNucleiBlob = im2bw(MSERNucleiBlob);
 
        % Form Blobs Set
        MSERNucleiBlobSet{i,1} = MSERNucleiBlob;
 
        t_mserblobs(i) = toc;
    end
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(outputpath, storageCommonPath, 'MSERNucleiBlobSet.mat'), ...
        'MSERNucleiBlobSet', 't_mserblobs');    
end

%% Using Rules to find nucleus from MSERNucleiBlobSet
%====================================================================
try
    load(strcat(outputpath, storageCommonPath, 'NucleiMask.mat'), ...
        'NucleiMaskSet');
    fprintf('Step 4 - Filter MSER Nuclei Blobs by Rules...done!\n');
catch
    fprintf('Step 4 - Filter MSER Nuclei Blobs by Rules...');
    NucleiMaskSet = cell(imNum,1);
    
    % Parameters for Rules
        maxNucleiEccentricity = 0.90;
        minNucleiArea = 100;
        maxNucleiClumpAreaRatio = 0.1;

    % timer
    t_mserclean = zeros(size(imNucleiSet,1),1);
    
    for i = 1:imNum
        tic;
        clumpMask = SceneCytoClumpMaskSet{i,1};
        MSERNucleiBlobMask = ~MSERNucleiBlobSet{i,1};
        nucleiMask_i = zeros(size(MSERNucleiBlobMask));
        
        MSERBlobStats = regionprops(MSERNucleiBlobMask, ...
            'Eccentricity', 'PixelIdxList', 'Area');
        for j = 1:length(MSERBlobStats)
            
            % Area Ratio of Nuclei and its Clump
            areaRatio = ComputeAreaRatio_NucleiClump(MSERBlobStats(j,1).Area, ...
                                                     MSERBlobStats(j,1).PixelIdxList,...
                                                     clumpMask);
            % Average Intensity of this MSER Blob
            meanIntensity = mean(imNucleiSet{i,1}(MSERBlobStats(j,1).PixelIdxList(:)));
            
            % Filter MSER Nuclei Blobs by Eccentricity and Area Ratio 
            if MSERBlobStats(j,1).Area > 300 && meanIntensity >= 200
                continue;
            end
            
            if MSERBlobStats(j,1).Eccentricity < maxNucleiEccentricity && ...
               MSERBlobStats(j,1).Area > minNucleiArea &&...
               areaRatio < maxNucleiClumpAreaRatio && ...
               meanIntensity <= 170
           
                nucleiMask_i(MSERBlobStats(j,1).PixelIdxList(:)) = 1;
                
            elseif MSERBlobStats(j,1).Area > 300 && ...
                   MSERBlobStats(j,1).Eccentricity < maxNucleiEccentricity && ...
                   areaRatio < maxNucleiClumpAreaRatio && ...
                   meanIntensity <= 200
                
                   nucleiMask_i(MSERBlobStats(j,1).PixelIdxList(:)) = 1;
            end
        end
        
        NucleiMaskSet{i,1} = nucleiMask_i;
        
        t_mserclean(i) = toc;
    end

    fprintf('\tStep 4 - DONE!\n');
    
    % Save section variable to file
    save(strcat(outputpath, storageCommonPath, 'NucleiMask.mat'),...
        'NucleiMaskSet', 't_mserclean');
end

