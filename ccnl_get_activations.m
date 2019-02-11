function activations = ccnl_get_activations(EXPT,model,mask,subjects)
    
    % Extract activation coefficients from a mask.
    % Caution: don't use for too many voxels
    %
    % USAGE: activations = ccnl_get_activations(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   activations - [nScans x nVoxels x nSubjects] activations
    %
    % Momchil Tomov, Aug 2018
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        %assert(isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');
       
        % extract voxel activations
        if strcmp(mask_format, 'cor')
            Y = spm_data_read(SPM.xY.VY, 'xyz', mask');
        else
            Y = spm_data_read(SPM.xY.VY, mask);
        end

        %-Whiten/Weight data and remove filter confounds
        KWY = spm_filter(SPM.xX.K,SPM.xX.W*Y);

        activations(:,:,s) = KWY;

        fprintf('Computed activations for subject %d\n', subj);
    end
