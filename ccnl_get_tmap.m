function tmap = ccnl_get_tmap(EXPT,model,regressor,mask,subjects)
    
    % Extract t-statistics from a mask. Also see ccnl_get_beta.
    % Works only for single-regressor contrasts.
    % You need to have run ccnl_fmri_con first.
    %
    % USAGE: tmap = ccnl_get_tmap(EXPT,model,regressor,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   regressor - regressor name or regressor index
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   tmap - [nSubjects x nVoxels] t-statistics
    %
    % Momchil Tomov, Aug 2018
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end
   
    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    if ischar(regressor)
        regressor = strtrim(regressor);
    end

    % load contrasts
    load(fullfile(EXPT.modeldir,['model',num2str(model)],'contrasts.mat'));

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);

        t = get_beta_or_tmap_helper(regressor, modeldir, Vmask, mask, contrasts, mask_format, 'spmT');
    
        if size(t, 1) > 1
            load(fullfile(modeldir,'SPM.mat'));
            b = get_beta_or_tmap_helper(regressor, modeldir, Vmask, mask, SPM.xX.name, mask_format, 'beta');
            assert(isequal(size(b), size(t)), 'Should have 1 t-map only or exactly 1 t-map per regressor');

            se = b ./ t;
            se = sqrt(nansum(se.^2)) / size(se, 1);
            b = nanmean(b, 1);
            t = b ./ se;
        end
        tmap(s,:) = t;
        
    end

end
