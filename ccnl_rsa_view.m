function ccnl_rsa_view(EXPT, rsa_idx, model_idx, T, all_subject_rhos, roi_masks)

    % View RSA group-level t-map from ccnl_rsa.
    % Requires bspmview.
    %
    % USAGE:
    %   ccnl_rsa_view(EXPT, rsa_idx, model_idx, T, roi_masks)
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   model_idx - which model to show 
    %   T - t-stats for each ROI, as returned by ccnl_rsa
    %   all_subject_rhos - rho's for all subjects, as returned by ccnl_rsa
    %   roi_masks -- ROI mask as passed to ccnl_rsa
    %
    % Momchil Tomov, May 2020 

    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);

    % create sample rsa
    rsa = EXPT.create_rsa(rsa_idx, 1);
    rsadir

    % initialize empty tmap
    V = spm_vol(rsa.mask);
    tmap = nan(V.dim);

    % hacks to make it save the t-map as a t-map
    V.fname = fullfile(rsadir, ['temp_tmap_', num2str(model_idx), '.nii']); % change immediately!
    V.dt = [16 0];
    V.private.dat.dtype = 'FLOAT32-LE';
    V.private.dat.fname = V.fname;

    df = size(all_subject_rhos, 3) - 1;
    for i = 1:length(roi_masks)
        roi_mask = roi_masks{i};

        % load roi mask
        [roi_mask_format, roi_mask] = get_mask_format_helper(roi_mask);
        assert(strcmp(roi_mask_format, 'mask'), 'Improper mask');

        tmap(roi_mask) = T(i, model_idx);
    end

    %V.descrip = sprintf('SPM{T_[%d.0]}', df); % hack to write the degrees of freedom, to allow thresholding in bspmview

    % save tmap
    V.fname
    spm_write_vol(V, tmap);

    % view tmap
    struc = fullfile(EXPT.modeldir,'mean.nii');
    if exist(struc,'file')
        bspmview(V.fname, struc);
    else
        bspmview(V.fname);
    end
end
