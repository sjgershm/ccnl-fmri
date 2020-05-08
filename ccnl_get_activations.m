function activations = ccnl_get_activations(EXPT,model,mask,subjects,apply_filter)
    
    % Extract activation coefficients from a mask.
    % Caution: don't use for too many voxels
    %
    % USAGE: activations = ccnl_get_activations(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number (any one will do)
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %   apply_filter (optional) - whether to high-pass filter data like SPM does (default false)
    %
    % OUTPUTS:
    %   activations - [nScans x nVoxels x nSubjects] activations
    %
    % Momchil Tomov, Aug 2018
   
    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('apply_filter', 'var')
        apply_filter = false;
    end

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

        if apply_filter
            % high-pass filter data
            Y = spm_filter(SPM.xX.K,Y);
        end

        activations(:,:,s) = Y;

        fprintf('Computed activations for subject %d\n', subj);
    end
