function residuals = ccnl_get_residuals(EXPT,model,mask,subjects)
    
    % Extract residuals coefficients from a mask.
    % Caution: don't use for too many voxels
    %
    % USAGE: residuals = ccnl_get_residuals(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in MNI coordinates as a [N x 3] matrix
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   residuals - [nScans x nVoxels x nSubjects] residuals
    %
    % Sam Gershman, Nov 2016
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if mask_format == 'mask' || islogical(mask)
        mask = find(mask);
    end

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
       
        % extract voxel activations
        nScans = length(SPM.xY.VY); % # of TRs
        for j=1:nScans
            if mask_format == 'mni'
                cor = mni2cor(mask, SPM.xY.VY(j).mat);
                Y(j,:) = spm_data_read(SPM.xY.VY(j), 'xyz', cor');
            else
                Y(j,:) = spm_data_read(SPM.xY.VY(j), mask);
            end
        end

        %-Whiten/Weight data and remove filter confounds
        KWY = spm_filter(SPM.xX.K,SPM.xX.W*Y);

        % extract residuals
        res = spm_sp('r',SPM.xX.xKXs,KWY);

        % sanity check
        ResSS = sum(res.^2); %-Residual SSQ
        V = spm_vol(fullfile(modeldir,'ResMS.nii'));
        if mask_format == 'mni'
            cor = mni2cor(mask, V.mat);
            ResMS = spm_data_read(V, 'xyz', cor');
        else
            ResMS = spm_data_read(V, mask);
        end
        assert(immse(ResSS / SPM.xX.trRV, ResMS) < 1e-9, ['Computed residuals don''t match ResMS.nii for subject', num2str(s)]);  % ResMS = ResSS scaled by tr(RV)

        residuals(:,:,s) = res;

        fprintf('Computed residuals for subject %d\n', s);
    end