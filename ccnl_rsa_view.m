function ccnl_rsa_view(EXPT, rsa_idx, model_idx)

    % View RSA group-level t-map from ccnl_rsa_searchlight.
    % Requires bspmview.
    %
    % USAGE:
    %   ccnl_rsa_view(EXPT, rsa_idx, model_idx)
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   model_idx - which model to show 
    %
    % Momchil Tomov, Oct 2018

    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);

    % create sample rsa
    rsa = EXPT.create_rsa(rsa_idx, 1);

    % initialize empty tmap
    V = spm_vol(rsa.mask);
    tmap = nan(V.dim);

    % hacks to make it save the t-map as a t-map
    V.fname = fullfile(rsadir, ['searchlight_tmap_', num2str(model_idx), '.nii']); % change immediately!
    V.dt = [16 0];
    V.private.dat.dtype = 'FLOAT32-LE';
    V.private.dat.fname = V.fname;

    % compute t-map, if it hasn't been computed already
    if ~exist(V.fname, 'file')
        df = NaN;
        files = dir(rsadir);
        for i = 1:length(files)
            file = files(i).name;
            if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
                disp(['Loading ', file]);
                load(fullfile(rsadir, file));

                for j = 1:size(cor, 1)
                    tmap(cor(j,1), cor(j,2), cor(j,3)) = T(j, model_idx);
                end
                df = size(all_subject_rhos, 3) - 1;
            end
        end

        V.descrip = sprintf('SPM{T_[%d.0]}', df); % hack to write the degrees of freedom, to allow thresholding in bspmview

        % save tmap
        V.fname
        spm_write_vol(V, tmap);
    end

    % view tmap
    struc = fullfile(EXPT.modeldir,'mean.nii');
    bspmview(V.fname, struc);
end
