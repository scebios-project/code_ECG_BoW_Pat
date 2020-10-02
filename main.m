clear all
close all

% Add local paths
addpath(genpath('encoding'));
addpath(genpath('pwmetric'));
addpath('LIBSVM\windows')
addpath(genpath('C:\Program Files\Matlab\R2018b\toolbox\WaveLab850'));
addpath(genpath('LCKSVD\OMPbox'));
addpath(genpath('LCKSVD\ksvdbox'));
addpath(genpath('LCKSVD'));
addpath(genpath('spams-matlab-v2.6\build'));
addpath(genpath('VLFeat\toolbox\mex'));

%======================
% Set parameters
%======================
basename  = 'ECG_2';
mkdir('data')
mkdir(['data/' basename])
mkdir(['data/' basename '/plot'])
params.Data_DIR     = 'data\ECG_data_immutable\';
params.Codebook_DIR = ['data\' basename '\'];
params.Feature_DIR  = ['data\' basename '\'];
params.Results_DIR  = ['data\' basename '\'];
params.Plot_DIR     = ['data\' basename '\plot\'];

% Prepare logfile for this run
starttime = datetime('now');
logfile = sprintf('%slog_%s.txt', params.Results_DIR, starttime);
logfile = strrep(logfile, ' ', '_');
logfile = strrep(logfile, ':', '_');
diary(logfile);
diary on;

fprintf('%s --- Script started\n',datetime('now'));

params.numClass=4;          % Number of classes
%
params.test_set_class_size  = 100;   % how many recordings to use for test, per class
params.train_set_class_size = 179;   % how many recordings to use for train, per class
%
params.sub_length  = 512;            % window length
params.inter_point = 32;             % window step
params.size_feat = 128;               % how many feature coeffs to keep
params.level = 7;                     % coarsest level of DWT 
%
%params.transfType='randproj';
%params.transfType='dct';
params.transfType='dwt';
%
params.codebookSize = 500;  
%params.codebookSize=400;
%params.codebookSize=40;
%
params.VQtype='sparse';     
params.dl_lambda = 0.5;
params.dl_iter   = 100;
%params.VQtype='LC-sparse';
%params.VQtype='sparseMultiD';
%params.VQtype='LC-KSVD';
%
params.encType = 'SPC';
params.poolType = 'sum-pooling';
%params.poolType = 'No';
%params.poolType = 'sum-pooling-size';
%params.poolType = 'max-pooling';
params.poolSize = 400;
params.poolStride = 200;
%params.normType = 'No';
params.normType = 'L2';                 % 'L2'; 'L1'; 'Power'; 'Power+L2'; 'Power+L1';
%params.normType = 'mean';                 % 'L2'; 'L1'; 'Power'; 'Power+L2'; 'Power+L1';
%
params.classDistType = [1 0 1];
%params.classType = 'NN';
params.classType = 'SVM';
%params.classType = 'LC-KSVD';
%params.sub_length = 128;            % window length
%params.inter_point = 4;             % window step
%


params.codebook_subsample_set_factor = 1;  % should we use less feature vectors for codebook training
%
params.global_iterations = 1;      % run everything this many times, for averaging the results
% Cleaning:
params.train_classes_to_clean = [];    % which classes to apply cleaning on
params.test_classes_to_clean  = [];    % which classes to apply cleaning on
%params.train_classes_to_clean = [1, 2, 3, 4];    % which classes to apply cleaning on
%params.test_classes_to_clean = [1, 2, 3, 4];    % which classes to apply cleaning on
%params.clean_type_detector = {'pan-tompkins', 'DWT', 'DWT', 'DWT'};
%params.clean_type_thresh = {'fixed', 'percent', 'percent', 'percent'};
%params.clean_dist_thresh = [0.5, 0.9, 0.9, 0.9];
%params.clean_type_detector = {'pan-tompkins', 'pan-tompkins', 'pan-tompkins', 'pan-tompkins'};
params.clean_type_detector = {'DWT', 'DWT', 'DWT', 'DWT'};
%params.clean_type_detector = {'pan-tompkins', 'pan-tompkins', 'pan-tompkins', 'pan-tompkins'};
params.clean_type_thresh = {'percent', 'percent', 'percent', 'percent'};
%params.clean_type_thresh = {'fixed', 'fixed', 'fixed', 'fixed'};
%params.clean_dist_thresh = [1, 1, 1, 1];
params.clean_dist_thresh = 1 * [1, 1, 1, 1];
params.clean_a_coeff  = 200;        
params.clean_a_offset = -100;
params.num_std_limit = 1;
%params.clean_function = @my_clean_ecg;
params.clean_function = @my_clean_ecg_new;
%
%params.tran_remove_mean = true;     % remove mean value from feature vectors, after feature extraction ?
params.tran_remove_mean = false;     % remove mean value from feature vectors, after feature extraction ?
%
params.plot = 0;

fprintf('%s --- Current parameters\n',datetime('now'));
disp(params);

% Initialize random number generator and save its state
rngstate_filename = [params.Feature_DIR, '0_rngstate.mat'];
if ~exist(rngstate_filename)
    rng('shuffle');
    rngstate = rng;
    fprintf('%s --- Randomly initialized RNG state, current state (JSON encoded):\n',datetime('now'));
    fprintf('%s\n',jsonencode(rngstate));
    fprintf('%s --- Saving RNG state to %s\n',datetime('now'), rngstate_filename);
    save(rngstate_filename, 'rngstate');
    changed_upstream = true;
else
    load(rngstate_filename);
    rng(rngstate);
    fprintf('%s --- Loaded RNG state from %s:\n',datetime('now'), rngstate_filename);
    fprintf('%s\n',jsonencode(rngstate));
end

if strcmp(params.transfType,'randproj')
    global Phi
    %Phi = randn(32,128);
    Phi = randn(params.size_feat, params.sub_length);
    %load([params.Data_DIR,'Phi.mat']);
end

% Set this to true in order to force rerun of everything from this point onwards
% (disable caching)
%changed_upstream = true;
changed_upstream = false;

%hit_rate_allexp = [];
conf_matrix_allexp = {};
F1_allexp = {};
for p = 1:params.global_iterations
    params.permIdx = p;
    fprintf('                      ----------\n');
    fprintf('%s Iteration %d / %d\n',datetime('now'), p, params.global_iterations);

    %================================================
    % Read data files and save data in matrix form
    %================================================
    data_save_filename = [params.Feature_DIR, '1_data_', num2str(p), '.mat'];
    if changed_upstream || ~exist(data_save_filename)
        fprintf('%s --- Reading raw data\n',datetime('now'));

        annotations = readtable([params.Data_DIR , '/training2017/REFERENCE.csv']);
        %annotations = readtable('data/ECG/training2017/REFERENCE-original.csv');

        names = annotations{:,1};
        cats = cell2mat(annotations{:,2});

        names_N = names(cats=='N');
        names_A = names(cats=='A');
        names_O = names(cats=='O');
        names_T = names(cats=='~');

        data_N = {};
        for i=1:numel(names_N)
            filename = [names_N{i} '.mat'];
            load([params.Data_DIR 'training2017/' filename]);
            data_N{i} = val';
        end

        data_A = {};
        for i=1:numel(names_A)
            filename = [names_A{i} '.mat'];
            load([params.Data_DIR 'training2017/' filename]);
            data_A{i} = val';
        end

        data_O = {};
        for i=1:numel(names_O)
            filename = [names_O{i} '.mat'];
            load([params.Data_DIR 'training2017/' filename]);
            data_O{i} = val';
        end

        data_T = {};
        for i=1:numel(names_T)
            filename = [names_T{i} '.mat'];
            load([params.Data_DIR 'training2017/' filename]);
            data_T{i} = val';
        end

        data  = {data_N,  data_A,  data_O,  data_T};
        names = {names_N, names_A, names_O, names_T};

        %======================================
        % Split test set and train set
        %======================================
        % Randomly split test set and train set
        test_set_class_size = params.test_set_class_size;  % test set size for a single class
        rp_all = {};
        for i = 1:numel(data)
            rp = randperm(size(data{i}, 2));
            rp_all{i} = rp;
            data_test{i}  =  data{i}(rp(1:test_set_class_size));
            names_test{i} = names{i}(rp(1:test_set_class_size));
            
            %data_train{i} = data{i}(rp(test_set_class_size+1 : test_set_class_size+train_set_class_size));
            data_train{i}  =  data{i}(rp(test_set_class_size+1 : end));
            names_train{i} = names{i}(rp(test_set_class_size+1 : end));
        end 

        fprintf('%s --- Saving raw data in %s\n',datetime('now'), data_save_filename);
        save(data_save_filename,'data', 'data_test', 'data_train', 'rp_all', 'names_test', 'names_train', '-v7.3');
        changed_upstream = true;  % recompute everything downstream
    else
        fprintf('%s --- Loading existing raw data from %s\n',datetime('now'), data_save_filename);
        load(data_save_filename,'data_test', 'data_train', 'rp_all', 'names_test', 'names_train');
    end

    [conf_matrix, F1]  = run_BoW(data_train, data_test, params, changed_upstream);
    %hit_rate_allexp = [hit_rate_allexp; hit_rate_all];
    conf_matrix_allexp = [conf_matrix_allexp; conf_matrix];
    F1_allexp = [conf_matrix_allexp; conf_matrix];
end

results_filename = sprintf('%sresults_%s.mat', params.Results_DIR, starttime);
results_filename = strrep(results_filename, ' ', '_');
results_filename = strrep(results_filename, ':', '_');

%fprintf('%s --- Results:\n',datetime('now'));
%disp(hit_rate_allexp);

%fprintf('%s --- Average:\n',datetime('now'));
%disp(mean(hit_rate_allexp));

fprintf('%s --- Saving classification results in %s\n',datetime('now'), results_filename);
save(results_filename, 'conf_matrix_allexp', 'F1_allexp', 'params', 'rngstate');

fprintf('%s - Finished\n',datetime('now'));
diary off