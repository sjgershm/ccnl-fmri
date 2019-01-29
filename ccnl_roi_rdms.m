function [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks, subjects)

    % Compute the neural RDMs for given ROIs. Also see ccnl_searchlight_rdms.m
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks [, subjects])
    %
    % EXAMPLE:
    %   [Neural] = ccnl_roi_rdms(exploration_expt(), 1, 'masks/hippocampus.nii')
    %   showRDMs(Neural(1).subj);
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   roi_masks - mask name or cell array of mask names. Could pass 3D masks instead.
    %   subjects - (optional) list of subjects
    %
    % OUTPUT:
    %   Neural - struct array, one element per sphere, with the fields:
    %      .name - mask name
    %      .event - within-trial event of the neural activation
    %      .n - number of voxels in sphere
    %      .subj - struct array with subject RDMs for given sphere; has following fields:
    %         .RDM - [nTrials x nTrials] RDM based on neural activation
    %   
    % Momchil Tomov, Sep 2018

    rng('shuffle'); % for parallelization

    if ~iscell(roi_masks)
        roi_masks = {roi_masks};
    end

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    % TODO dedupe with ccnl_searchlight_rdms.m

    % create rsa folder if none exists
    if ~isdir(EXPT.rsadir); mkdir(EXPT.rsadir); end
    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);
    if ~isdir(rsadir); mkdir(rsadir); end

    % load example rsa
    rsa = EXPT.create_rsa(rsa_idx, 1);
    mask = rsa.mask;

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');

    % for each subject
    for s = 1:length(subjects)
        subj = subjects(s);

        rsa = EXPT.create_rsa(rsa_idx, subj);

        betas_filename = fullfile(rsadir, sprintf('betas_%d.mat', subj));
        disp(betas_filename);

        % load (cached) betas
        if ~exist(betas_filename, 'file')
            tic
            disp('loading betas from .nii files...');

            % load betas 
            modeldir = fullfile(EXPT.modeldir,['model',num2str(rsa.glmodel)],['subj',num2str(subj)]);
            load(fullfile(modeldir,'SPM.mat'));
            which = contains(SPM.xX.name, rsa.event); % betas for given event
            which(which) = rsa.which_betas; % of those, only betas for given trials
            cdir = pwd;
            cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
            B = spm_data_read(SPM.Vbeta(which), find(mask));
            cd(cdir);

            % save file in "lock-free" fashion
            % b/c parallel jobs might be doing the same
            tmp_filename = [betas_filename, random_string()];
            save(tmp_filename, 'B', '-v7.3');
            movefile(tmp_filename, betas_filename); % TODO assumes this is instantaneous

            toc
        else
            tic
            disp('loading cached betas from .mat file...');
            load(betas_filename);
            toc
        end

        tic
        disp('computing RDMs...');

        % for each voxel
        for i = 1:length(roi_masks)
            roi_mask = roi_masks{i};

            if ischar(roi_mask)
                mask_name = roi_mask;
            else
                mask_name = sprintf('mask_%d', i);
            end
            disp(mask_name);

            % load roi mask
            [roi_mask_format, roi_mask] = get_mask_format_helper(roi_mask);
            assert(strcmp(roi_mask_format, 'mask'), 'Improper mask');

            roi_mask = roi_mask & mask; % exclude out-of-brain voxels
            roi_mask(mask) = roi_mask(mask) & ~any(isnan(B), 1)'; % exclude nan voxels

            % compute RDM
            if sum(roi_mask(:)) == 0
                % sometimes (rarely) they're all NaNs
                Neural(i).subj(s).RDM = [];
            else
                Neural(i).subj(s).RDM = squareRDMs(pdist(B(:, find(roi_mask(mask))), 'cosine'));
            end
            assert(sum(any(isnan(Neural(i).subj(s).RDM))) == 0, 'Found NaNs in RDM -- should never happen');

            % metadata
            Neural(i).subj(s).id = subj;
            Neural(i).name = mask_name;
            Neural(i).subj(s).name = Neural(i).name;
            Neural(i).event = rsa.event;
            Neural(i).num_voxels = sum(roi_mask(:));
        end

        toc
    end

end

