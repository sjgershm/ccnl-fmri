function beta = ccnl_get_beta(EXPT,model,regressor,mask,subjects)
    
    % Extract beta coefficients from a mask.
    %
    % USAGE: beta = ccnl_get_beta(EXPT,model,regressor,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   regressor - regressor name or regressor index
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in MNI coordinates as a [N x 3] matrix
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   beta - [nSubjects x nVoxels] beta coefficients
    %
    % Sam Gershman, Nov 2016
    % Momchil Tomov, Aug 2018
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    if ischar(regressor)
        regressor = strtrim(regressor);
    end

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
       
        b = get_beta_or_tmap_helper(regressor, modeldir, Vmask, mask, SPM.xX.name, mask_format, 'beta');
        
        beta(s,:) = nanmean(b, 1);
        
    end
