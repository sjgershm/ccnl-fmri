function bic = ccnl_bic(EXPT,model,mask,subjects)
    
    % Compute Bayesian information criterion (BIC) for a model.
    %
    % USAGE: bic = ccnl_bic(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   mask - a mask image name (e.g., 'mask.nii'), a set of voxel
    %          indices, or a binary vector
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   bic - [nSubjects x 1] vector of BIC values
    %
    % Sam Gershman, Oct 2016
    
    if nargin < 3; subjects = 1:length(EXPT.subject); end
    
    % load mask
    if ischar(mask)
        V = spm_vol(mask);
        Y = spm_read_vols(V);
        mask = Y~=0 & ~isnan(Y);
    end
    
    bic = zeros(length(subjects),1);
    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        [N,K] = size(SPM.xX.X); % number of datapoints and parameters
        V = spm_vol(fullfile(modeldir,'ResMS.nii'));    % residual variance image
        Y = spm_read_vols(V);
        bic(s) = N*sum(log(Y(mask))) + K*log(N);
    end