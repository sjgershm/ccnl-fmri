function [Rho, H, T, P, all_subject_rhos, Behavioral, Neural] = ccnl_rsa(EXPT, rsa_idx, roi_masks, subjects)

    % RSA for given ROIs. Also see ccnl_rsa_searchlight.m
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [Rho, H, T, P, all_subject_rhos, Behavioral, Neural] = ccnl_rsa(EXPT, rsa_idx, roi_masks)    
    %
    % EXAMPLE:
    %   [Rho, ~, T, P] = ccnl_rsa(exploration_expt(), 1, 'masks/hippocampus.nii');
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   roi_masks - mask name or cell array of mask names. Could pass 3D masks instead.
    %   subjects - (optional) list of subjects
    %
    % OUTPUT:
    %   Rho - [nROIs x nModels] matrix of Spearman rank correlations (averaged across subjects)
    %   H - [nROIs x nModels] matrix of hypothesis outcomes
    %   T - [nROIs x nModels] matrix of t-statistics
    %   P - [nROIs x nModels] matrix of p-values
    %   all_subject_rhos - [nROIs x nModels x nSubjects] matrix of Spearman rank correlations
    %   Behavioral - struct array with behavioral RDMs (see ccnl_behavioral_rdms.m)
    %   Neural - struct array with neural RDMs (see ccnl_roi_rdms.m)
    %
    % Momchil Tomov, Oct 2018

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    % create rsa folder if none exists
    if ~isdir(EXPT.rsadir); mkdir(EXPT.rsadir); end
    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);
    if ~isdir(rsadir); mkdir(rsadir); end

    % get behavioral (model) RDMs
    [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx, subjects);

    % get searchlight RDMs
    [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks, subjects);

    % compute second-order correlations (similarity match)
    [Rho, H, T, P, all_subject_rhos] = ccnl_match_rdms(Neural, Behavioral, control);

