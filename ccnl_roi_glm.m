function [beta, ResMS] = ccnl_roi_glm(EXPT, model, mask, regressor, subjects)

    % Run a GLM using average activity within a ROI. Uses logic from SPM12.
    % You need to have run ccnl_fmri_glm first. 
    % Caution: don't use for too many voxels
    %
    % USAGE:
    %   [beta, ResMS] = ccnl_roi_glm(EXPT, model, mask, regressor, subjects)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   regressor - regressor name or regressor index
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          or a binary vector/mask in native space.
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUT:
    %   beta - [nSubjects] beta coefficients
    %   ResMS - [nSubjects] residual variance
    %
    % Momchil Tomov, Oct 2018

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end
    
    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    regressor = strtrim(regressor);

    for s = 1:length(subjects)
        fprintf('Subject %d\n', s);

        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');

        nScans = length(SPM.xY.VY); % # of TRs

        % extract voxel activations
        Y = spm_data_read(SPM.xY.VY, mask);
        
        % average activation within ROI
        Y = mean(Y,2);

        %-Whiten/Weight data and remove filter confounds
        % TODO look into this
        KWY = spm_filter(SPM.xX.K,SPM.xX.W*Y);

        % extract residuals
        res = spm_sp('r',SPM.xX.xKXs,KWY);
        ResSS = sum(res.^2); %-Residual SSQ
        
        %-Weighted Least Squares estimation
        %======================================================================
        b         = SPM.xX.pKX*KWY;                     %-Parameter estimates

        % get the beta we care about
        idx = [];
        for i = 1:length(SPM.xX.name)
            if ~isempty(strfind(SPM.xX.name{i},[regressor,'*'])) || ~isempty(strfind(SPM.xX.name{i},[regressor,'^'])) || endsWith(SPM.xX.name{i},regressor)
                idx = [idx i];
            end
        end

        beta(s) = mean(b(idx));
        ResMS(s) = ResSS / SPM.xX.trRV;
    end
