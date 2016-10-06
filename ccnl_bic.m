function bic = ccnl_bic(EXPT,model,mask,subjects)
    
    % Compute Bayesian information criterion (BIC) for a model.
    %
    % USAGE: bic = ccnl_bic(EXPT,model,mask,[subjects])
    
    if nargin < 3; subjects = 1:length(EXPT.subject); end
    
    % load mask
    V = spm_vol(mask);
    Y = spm_read_vols;
    ix = Y~=0 && ~isnan(Y);
    
    bic = zeros(length(subjects),1);
    for s = 1:length(subjects)
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        [N,K] = size(SPM.xX.X); % number of datapoints and parameters
        V = spm_vol(fullfile(modeldir,'ResMS.nii'));    % residual variance image
        Y = spm_read_vols(V);
        bic(s) = N*sum(log(Y(ix))) + K*log(N);
    end