function [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj, substring, mask)

    % Extracts beta coefficients for all regressors that contain a given substring.
    % Useful for MVPA, beta series correlation and PPI.
    % Note that it's not recommended if events are closer than ~8 s together (Poldrack's 2011 handbook of fMRI, section 8.2.3)
    %
    % USAGE: 
    %   B = ccnl_get_beta_series(EXPT, model, subj, substring, mask)
    %   b = mean(B, 1); % average across trials / runs -- useful for e.g. inter-subject RSA,
    %   % ...or...
    %   b = mean(B, 2); % average across voxels -- useful for e.g. PPI or beta series correlation.
    %
    % EXAMPLE:
    %   VS = ccnl_get_beta_series(exploration_expt(), 57, 1, 'trial_onset', 'masks/NAC.nii');
    %   VS = nanmean(VS, 2);
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subject - which subject to analyze (one at a time because # of trials might differ across subjects)
    %   substring - regressor substring, e.g. 'trial_onset_'
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %
    % OUTPUTS:
    %   beta - [nTrials x nVoxels] beta coefficients
    %   names - {nTrials x 1} cell array of regressor names
    %
    % Momchil Tomov, Apr 2019

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');

    % load betas 
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    assert(isempty(Vmask) || isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');
    assert(ndims(mask) < 3 || isequal(SPM.Vbeta(1).dim, size(mask)), 'Different dimensions between mask and betas');

    which = contains(SPM.xX.name, substring); % betas for given event
    cdir = pwd;
    cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
    B = spm_data_read(SPM.Vbeta(which), find(mask));
    cd(cdir);

    names = SPM.xX.name(which)';
           
    assert(size(B,1) > 0, 'no betas - likely wrong substring');
