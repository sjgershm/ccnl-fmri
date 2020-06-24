function [Neural, cor] = ccnl_searchlight_rdms(EXPT, rsa_idx, inds, radius, subjects, distance_metric)

    % Compute the searchlight neural RDMs
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [Neural] = ccnl_searchlight_rdms(EXPT, rsa_idx, inds)
    %
    % EXAMPLE:
    %   [Neural] = ccnl_searchlight_rdms(exploration_expt(), 1, 1:100)
    %   showRDMs(Neural(1).subj);
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   inds - which voxels to use for sphere centers, as indices in mask(:)
    %   radius - (optional) searchlight radius in voxels
    %   subjects - (optional) list of subjects
    %
    % OUTPUT:
    %   Neural - struct array, one element per sphere, with the fields:
    %      .name - sphere name
    %      .cor - coordinates of sphere center
    %      .mni - MNI coordinates of sphere center
    %      .radius - sphere radius
    %      .event - within-trial event of the neural activation
    %      .n - number of voxels in sphere
    %      .subj - struct array with subject RDMs for given sphere; has following fields:
    %         .RDM - [nTrials x nTrials] RDM based on neural activation
    %   
    % Momchil Tomov, Sep 2018

    rng('shuffle'); % for parallelization

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('distance_metric', 'var') 
        distance_metric = 'cosine';
    end

    % create rsa folder if none exists
    if ~isdir(EXPT.rsadir); mkdir(EXPT.rsadir); end
    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);
    if ~isdir(rsadir); mkdir(rsadir); end

    % load example rsa
    rsa = EXPT.create_rsa(rsa_idx, 1);
    mask = rsa.mask;

    if ~exist('radius', 'var')
        radius = rsa.radius;
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');

    % get all voxel coordinates
    [x,y,z] = ind2sub(size(mask), find(mask));
    min_x = min(x);
    max_x = max(x);
    min_y = min(y);
    max_y = max(y);
    min_z = min(z);
    max_z = max(z);

    % get voxels we care about
    inds = inds(inds <= length(x));
    x = x(inds);
    y = y(inds);
    z = z(inds);
    cor = [x y z];

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
            disp('loading betas from .nii files...');

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
        for i = 1:length(x)
            % create searchlight mask
            sphere_mask = create_spherical_mask_helper(mask, x(i), y(i), z(i), radius, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);
            sphere_mask = sphere_mask & mask; % exclude out-of-brain voxels
            sphere_mask(mask) = sphere_mask(mask) & ~any(isnan(B), 1)'; % exclude nan voxels

            % compute RDM
            if sum(sphere_mask(:)) == 0
                % sometimes (rarely) they're all NaNs
                Neural(i).subj(s).RDM = [];
            else
                Neural(i).subj(s).RDM = squareRDMs(pdist(B(:, find(sphere_mask(mask))), distance_metric));
            end
            assert(sum(any(isnan(Neural(i).subj(s).RDM))) == 0, 'Found NaNs in RDM -- should never happen');
          
            % metadata
            Neural(i).subj(s).id = subj;
            Neural(i).cor = [x(i) y(i) z(i)];
            Neural(i).mni = cor2mni(Neural(i).cor, Vmask.mat);
            Neural(i).name = ['sphere_', sprintf('%d_%d_%d', Neural(i).mni), '_', rsa.event];
            Neural(i).subj(s).name = Neural(i).name;
            Neural(i).radius = radius;
            %Neural(i).event = rsa.event;
            Neural(i).idx = inds(i);
            Neural(i).num_voxels = sum(sphere_mask(:));
        end

        toc
    end

end


function [sphere_mask, sphere_coords] = create_spherical_mask_helper(mask, x, y, z, r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask)

    % does the bulk of the work in create_spherical_mask
    % convenient to use if getting multiple masks in bulk (saves some overhead)
    % Note that the coordinates are in native space, NOT in MNI space!
    %
    sphere_coords = [];

    sphere_mask = zeros(size(mask));

    for newx = floor(x - r) : ceil(x + r)
        if newx < min_x || newx > max_x, continue; end
        for newy = floor(y - r) : ceil(y + r)
            if newy < min_y || newy > max_y, continue; end
            for newz = floor(z - r) : ceil(z + r)
                if newz < min_z || newz > max_z, continue; end
                if ~mask(newx, newy, newz), continue; end
                if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > r^2, continue; end
                sphere_mask(newx, newy, newz) = 1;
               % mni = cor2mni([newx newy newz], Vmask.mat);
               % sphere_coords = [sphere_coords; mni];
                sphere_coords = [sphere_coords; newx newy newz];
            end
        end
    end

    sphere_mask = logical(sphere_mask);

end
