% SCRIPT FILE: 
%
%% 
tidy; 
warning('off', 'all');
[dn,ds] = loadnames('macros', chooseplatform);

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
