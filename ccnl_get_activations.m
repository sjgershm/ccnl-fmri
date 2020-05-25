function [activations, runs] = ccnl_get_activations(EXPT,model,mask,subjects,whiten,filter)
    
    % Extract activation coefficients from a mask.
    % Caution: don't use for too many voxels
    %
    % USAGE: [activations, runs] = ccnl_get_activations(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number (matters only if whitening, otherwise any one will do)
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %   whiten (optional) - whether to whiten data like SPM does (default true, like SPM)
    %   filter (optional) - whether to high-pass filter data like SPM does (default true, like SPM)
    %
    % OUTPUTS:
    %   activations{s} - [nScans x nVoxels] activations for subject s
    %   runs{s} - [nScans x 1] run/session IDs for subject s
    %
    % Momchil Tomov, Aug 2018
   
    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('whiten', 'var')
        whiten = true;
    end

    if ~exist('filter', 'var')
        filter = true;
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
        assert(isempty(Vmask) || isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');
       
        % extract voxel activations
        if strcmp(mask_format, 'cor')
            Y = spm_data_read(SPM.xY.VY, 'xyz', mask');
        else
            Y = spm_data_read(SPM.xY.VY, mask);
        end

        % see spm_spm.m
        if whiten
            % whiten data
            Y = SPM.xX.W*Y;
        end
        if filter
            % high-pass filter data
            Y = spm_filter(SPM.xX.K,Y);
        end

        activations{s} = Y;

        for r = 1:length(SPM.Sess)
            runs{s}(SPM.Sess(r).row,:) = r;
        end

        fprintf('Computed activations for subject %d\n', subj);
    end
