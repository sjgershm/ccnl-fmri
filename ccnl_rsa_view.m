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
    [~, V, tmap] = load_mask(rsa.mask);
    tmap(:) = NaN; % clear
    V.fname = fullfile(rsadir, ['searchlight_tmap_', num2str(model_idx), '.nii']); % change immediately!

    % compute t-map, if it hasn't been computed already
    if ~exist(V.fname, 'file')
        files = dir(rsadir);
        for i = 1:length(files)
            file = files(i).name;
            if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
                disp(['Loading ', file]);
                load(fullfile(rsadir, file));

                for j = 1:size(cor, 1)
                    tmap(cor(j,1), cor(j,2), cor(j,3)) = T(j, model_idx);
                end
            end
        end

        % save tmap
        spm_write_vol(V, tmap);
    end

    % view tmap
    struc = fullfile(EXPT.modeldir,'mean.nii');
    bspmview(V.fname, struc);
end
