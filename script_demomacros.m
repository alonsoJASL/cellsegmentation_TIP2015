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
    [X, xatt] = readParseSome(strcat(dn,ds),names);
    imNum = size(X,4);
else 
    [X, xatt] = readParseInput(strcat(dn,ds, 'man00', num2str(names),'.tif'));
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
                        
%% 

load(strcat(outputpath,storageCommonPath,'SceneCytoClumpMaskSet.mat'));
load(strcat(outputpath,storageCommonPath,'NucleiMask.mat'));
