function v = ccnl_vbm(EXPT, mask, subjects, tissue, modulated)

    % Voxel-based morphometry (VBM).
    % Extract (expected) number of voxels of given tissue type in ROI.
    %
    % USAGE: ccnl_vbm(EXPT, mask, [subjects], [tissue], [modulated])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - vector of subject numbers (default: all subjects)
    %   tissue (optional) - 'grey', 'white', or 'csf' (defaults to 'grey')
    %   modulated (optional) - whether to account for volume expansion/shrinkage during spatial normalization (defaults to true)
    %
    % OUTPUT:
    %    v - [nSubjects] number of voxels for each subject
    %
    % Momchil Tomov, Nov 2018

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end
    if ~exist('tissue', 'var')
        tissue = 'grey';
    end
    if ~exist('modulated', 'var')
        modulated = true;
    end

    if modulated
        prefix = 'mwc';
    else
        prefix = 'wc';
    end

    switch tissue
        case 'grey'
            prefix = [prefix, '1'];
        case 'white'
            prefix = [prefix, '2'];
        case 'csf'
            prefix = [prefix, '3'];
        otherwise
            assert(false, 'tissue must be grey, white, or csf');
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    % Loop over subjects
    for s = 1:length(subjects)
        subj = subjects(s);
        
        S = EXPT.subject(subj); % subject structure
        tissuefile = fullfile(S.datadir, [prefix, S.structural]);

        V = spm_vol(tissuefile);
        assert(isequal(V.dim, Vmask.dim), 'Different dimensions between mask and tissue volume');

        % extract tissue map
        if strcmp(mask_format, 'cor')
            Y = spm_data_read(V, 'xyz', mask');
        else
            Y = spm_data_read(V, mask);
        end

        v(s) = nansum(Y);
    end
        

end
