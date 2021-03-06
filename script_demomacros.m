% Script file: demo - Macrophages dataset
%
% 

%% Initialisation Linux

clear all;
close all;
clc;
warning('off', 'all');

plat = 'linux';
[dn,ds] = loadnames('macros', plat);
names = [55];% 75 95];

name = 'man00';

outputpath = fullfile(dn, [ds(1:end-1) '_levelset_OUT']);
storageCommonPath = 'Common/';
storageInitial= 'Initial/';
storageDist= 'LSF/';

if ~isdir(outputpath)
    mkdir(outputpath);
    mkdir(fullfile(outputpath, storageCommonPath));
    mkdir(fullfile(outputpath, storageInitial));
    mkdir(fullfile(outputpath, storageDist));
end

addpath(genpath(pwd));
%% Initialisation OSX 

clear all;
close all;
clc;

warning('off', 'all');

plat = 'osx';
[dn,ds] = loadnames('macros', plat);
names = [55 75 95];

name = 'man00';

outputpath = strcat(dn, ds(1:end-1), '_levelset_OUT/');
storageCommonPath = 'Common/';
storageInitial= 'Initial/';
storageDist= 'LSF/';
storageVoronoi = 'Voronoi/';

if ~isdir(outputpath)
    mkdir(outputpath);
    mkdir(strcat(outputpath, storageCommonPath));
    mkdir(strcat(outputpath, storageInitial));
    mkdir(strcat(outputpath, storageDist));    
    mkdir(strcat(outputpath, storageVoronoi));
end

addpath(genpath(pwd));

%% Initialisation WINDOWS 

clear all;
close all;
clc;

warning('off', 'all');

plat = 'win';
[dn,ds] = loadnames('macros', plat);
names = 55;%[55 75 95];

name = 'man00';

outputpath = strcat(dn, ds(1:end-1), '_levelset_OUT\');
storageCommonPath = 'Common\';
storageInitial= 'Initial\';
storageDist= 'LSF\';
storageVoronoi = 'Voronoi\';

if ~isdir(outputpath)
    mkdir(outputpath);
    mkdir(strcat(outputpath, storageCommonPath));
    mkdir(strcat(outputpath, storageInitial));
    mkdir(strcat(outputpath, storageDist));    
    mkdir(strcat(outputpath, storageVoronoi));
end

addpath(genpath(pwd));
%% Clear previous results and outputs.
% ======================================
unix(['rm -rf ' 32 outputpath]);

if ~isdir(outputpath)
    mkdir(outputpath);
    mkdir(strcat(outputpath, storageCommonPath));
    mkdir(strcat(outputpath, storageInitial));
    mkdir(strcat(outputpath, storageDist));
    mkdir(strcat(outputpath, storageVoronoi));
end

%% Parameters Initialisation for methods
%=======================================

% Distance Map Parameter 
beta_logistic_set = (5)';
% Joint Level Set Parameter
kappa_set = (13)';
chi_set = (3)';
% Joint Level Set Updating Policy (create necessary 
loop = 5;
if ~isdir(strcat(outputpath, storageDist, 'LSF1/'))
    for i = 1:loop
        mkdir(strcat(outputpath, storageDist, 'LSF', num2str(i), '/'));
    end
end

%========================================================
% Common parameters for GMM-based clump training/testing
inParam.minCellSize = 1000;
 
% Level set paramters for clump refinement (after GMM)
inParam.iter_in = 2;
inParam.iter_out = 3;
inParam.alfaGMM = -2.5;
inParam.lambdaGMM = 5;
 
% Level set parameters for nuclei refinement (after GMM)
inParam.iter_in_rawNuclei = 10;
inParam.iter_out_rawNuclei = 2;
inParam.alfa_rawNuclei = 5;
inParam.lambda_rawNuclei = 4;

%   Level set for segmentation
inParam.iter_in_extent = 20;
inParam.iter_out_extent = 2;

%     Radiating distance map
% nLinePartition = 10;
inParam.min_PtDistance = 5;
inParam.min_TinyFragments_DistMap = 500;

%% PURE GREYSCALE - Loading the images into memory

if length(names) > 1
    [X, xatt] = readParseSome(strcat(dn,ds),names);
    imNum = size(X,4);
else 
    [X, xatt] = readParseInput(strcat(dn,ds, 'man00', num2str(names),'.tif'));
    imNum = 1;
end

% create greyscale images (downscaled)
imCytoSet = cell(imNum,1);
for i=1:imNum
    Iaux = rgb2gray(X(:,:,:,i));
    Iaux = abs(Iaux - max(Iaux(:)));
    
    % for this method, a uint8 image is needed.
    % Rescaling is optional.
    imCytoSet{i} = im2uint8(imresize(Iaux,1));
end

imNucleiSet = imCytoSet;

clear Imaux;

%% RED AND GREEN SEPARATE - Loading the images into memory

if length(names) > 1
    [X, xatt] = readParseSome(filename,names);
    imNum = size(X,4);
else 
    filename = strcat('man00', num2str(names),'.tif');
    [X, xatt] = readParseInput(strcat(dn,ds,filename));
    imNum = 1;
end

% create greyscale images (downscaled)
imCytoSet = cell(imNum,1);
for i=1:imNum
    IRaux = X(:,:,1,i);
    IRaux = abs(IRaux - max(IRaux(:)));
    
    IGaux = rgb2gray(X(:,:,:,i));
    IGaux = abs(IGaux - max(IGaux(:)));
    
    % for this method, a uint8 image is needed.
    % Rescaling is optional.
    imNucleiSet{i} = im2uint8(imresize(IRaux,1));
    imCytoSet{i} = im2uint8(imresize(IGaux,1));
end

imNum = size(imCytoSet,1);

clear IGaux IRaux;

%% RUN STANDARD SEGMENTATION - Only find clumps and nuclei blobs
%
clc;
tic;
fullRunInSteps(imCytoSet, imNucleiSet, storageCommonPath, ...
                inParam,outputpath);
t = toc;
fprintf('\n FULL TIME CLUMP AND NUCLEI SEGMENTATION %5.3f.\n',t);

clear t;

%% RUN VORONOI SEGMENTATION - Low tech, but reliable!

tic; 
dataLa = cell(imNum,1);
for i=1:imNum
    [dataLa{i}, labAtt] = voronoiSegmentation(X(:,:,:,i), xatt(i));
end
toc;

% Need to store SceneCytoClumpMaskSet.mat and NucleiMask.m from dataLa{:}

NucleiMaskSet = cell(imNum,1);
SceneCytoClumpMaskSet = cell(imNum,1);

for i=1:imNum 
    NucleiMaskSet{i} = dataLa{i}(:,:,1)>0;
    SceneCytoClumpMaskSet{i} = dataLa{i}(:,:,2)>0;
end

save(strcat(outputpath,storageCommonPath, 'NucleiMask.mat'),'NucleiMaskSet');
save(strcat(outputpath,storageCommonPath, 'SceneCytoClumpMaskSet.mat'),...
    'SceneCytoClumpMaskSet');
save(strcat(outputpath, storageVoronoi, 'CompleteSegmentation.mat'), 'dataLa');


%% RUN LEVEL SET METHOD FOR OVERLAPPING
clc;

fullOverlappingSegmentation(imCytoSet, storageCommonPath, ...
                        inParam, outputpath, storageInitial, storageDist,...
                            beta_logistic_set, kappa_set, chi_set, loop);
                        
%% COMPARE RESULTS - Single image experiment

% Load first result form run.
close all

load(strcat(outputpath,storageDist, 'LSF5\',...
    'LSF_5_beta_5_kappa_13_chi_3_iterIn_20_iterOut_2.mat'));

resultlsf = LSF_5{1};
numCells = size(resultlsf,1);

height = size(resultlsf{1,1},1);
width = size(resultlsf{1,1},2);

overlappedCells = zeros(height, width, numCells);

n = 1;
for ix=1:numCells
    n = getPrimes(n,true);
    overlappedCells(:,:,ix) = n.*resultlsf{ix,1};
    
    %imagesc(overlappedCells(:,:,ix)); title(num2str(n));
    %pause;
end 



figure(21)
load(strcat(outputpath,storageCommonPath,'SceneCytoClumpMaskSet.mat'));
load(strcat(outputpath,storageCommonPath,'NucleiMask.mat'));
imagesc(SceneCytoClumpMaskSet{1}+NucleiMaskSet{1});
colormap bone;
title('LSM initial clumps and nuclei')

% Voronoi result
load(strcat(outputpath,storageVoronoi,'CompleteSegmentation.mat'));
% load ground truth
load(strcat(dn,ds(1:end-1),'_GT/', filename(1:end-4),'.mat'));

figure(1)
subplot(131);
imagesc(dataBin(:,:,2));
cooljet3;
title('Ground truth');

%figure(2)
subplot(132);
imagesc(changeOverlapRepresentation(overlappedCells)); 
cooljet3;
title('Segmentation Result');

%figure(3)
subplot(133);
imagesc(dataLa{1}(:,:,2));
cooljet3;
title('Voronoi');

% just to check 
figure(4)
imshow(X);

% initial guess 
load(strcat(outputpath, storageInitial,'initialPhis_beta_5.mat'));
initialGuess = roughInitialGuessMaskSet_byImage{1};

figure(22);
layer = 20;
imagesc(initialGuess{layer,1});
title(strcat('Initial guess for Layer: ', num2str(layer)));

figure(23)
imagesc(overlappedCells(:,:,layer));
title(strcat('Layer: ',num2str(layer))); 
%% COMPARISON AGAINST GROUND TRUTH 
% Voronoi vs. LSM

Blsm = cat(3, NucleiMaskSet{1}, changeOverlapRepresentation(overlappedCells),...
    zeros(size(NucleiMaskSet{1})));

[jmLSM, compattLSM] = segmentationComparison(Blsm,dataBin);
[bstats, gtstatsLSM] = getSegmentationStats(jmLSM, compattLSM);

sGT = compattLSM.stackedGT;
for i=1:size(sGT,3)
    sGT(:,:,i) = (sGT(:,:,i)>0).*gtstatsLSM.precision(i);
end

Bvoronoi = dataLa{1};

[jmV, compattV] = segmentationComparison(Bvoronoi,dataBin);
[bstatsV, gtstatsV] = getSegmentationStats(jmV, compattV);

sGT = compattV.stackedGT;
for i=1:size(sGT,3)
    sGT(:,:,i) = (sGT(:,:,i)>0).*gtstatsV.precision(i);
end

%%
figure(901)
subplot(121);
imagesc(max(sGT,[],3));
title('Segmentation precision on Level Set Method');
axis off;
colormap(1-hot);
colorbar;

%figure(902)
subplot(122);
imagesc(max(sGT,[],3));
title('Segmentation precision on V');
axis off;
colormap(1-hot);
colorbar;

%% 

load(strcat(outputpath,storageCommonPath,'SceneCytoClumpMaskSet.mat'));
load(strcat(outputpath,storageCommonPath,'NucleiMask.mat'));
