function [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks, subjects, return_B, distance_metric)

    % Compute the neural RDMs for given ROIs. Also see ccnl_searchlight_rdms.m
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks [, subjects], [return_B], [distance_metric])
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
    %   return_B - (optional) whether to return neural activity (BOLD or betas) used to generate RDM as part of Neural (default to false because it's too big)
    %   distance_metric - (optional) neural distance metric (default: cosine)
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

    if ~exist('distance_metric', 'var') 
        distance_metric = 'cosine';
    end

    if ~exist('return_B', 'var')
        return_B = false;
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

        if rsa.use_beta_series
            filename = fullfile(rsadir, sprintf('betas_%d.mat', subj));
        else
            filename = fullfile(rsadir, sprintf('BOLD_%d.mat', subj));
        end
        disp(filename);

        % load (cached) betas
        if ~exist(filename, 'file')
            tic
            disp('loading betas or BOLD from .nii files...');

            if rsa.use_beta_series
                % load betas
                %
                B = ccnl_get_beta_series(EXPT, rsa.glmodel, subj, rsa.event, rsa.mask);
            else
                % load BOLD
                %
                B = ccnl_get_activations(EXPT, rsa.glmodel, rsa.mask, subj, true, true); % whiten & filter; see Diedrichsen et al. 2016
                B = B{1};
            end

            % save file in "lock-free" fashion
            % b/c parallel jobs might be doing the same
            tmp_filename = [filename, random_string()];
            save(tmp_filename, 'B', '-v7.3');
            movefile(tmp_filename, filename); % TODO assumes this is instantaneous

            toc
        else
            tic
            disp('loading cached betas from .mat file...');
            load(filename);
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
                Neural(i).subj(s).RDM = squareRDMs(pdist(B(:, find(roi_mask(mask))), distance_metric));
            end
            assert(sum(any(isnan(Neural(i).subj(s).RDM))) == 0, 'Found NaNs in RDM -- should never happen');

            % metadata
            Neural(i).subj(s).id = subj;
            Neural(i).name = mask_name;
            Neural(i).subj(s).name = Neural(i).name;
            %Neural(i).event = rsa.event;
            Neural(i).num_voxels = sum(roi_mask(:));
            if return_B
                Neural(i).subj(s).B = B(:, find(roi_mask(mask)));
            end
        end

        toc
    end

end

