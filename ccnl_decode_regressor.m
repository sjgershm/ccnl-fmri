function dec = ccnl_decode_regressor(EXPT, glmodel, regressor, mask, lambda, subjects)

    % Invert the GLM (using ridge regression) to decode a regressor from the neural data.
    %
    % The formula is:
    %    dec = (y - sum x_i * b_i) * b_reg / (b_reg^2 + lambda)
    % where:
    %    dec = decoded regressor
    %    y = neural data
    %    x_i = regressor i
    %    b_i = beta coefficient for regressor i
    %    the sum is over all regressors i that are not our regressor
    %    b_reg = is the beta coefficient for our regressor
    %    lambda = ridge regularization parameter
    %
    % USAGE:
    %   dec = ccnl_decode_regressor(EXPT, glmodel, regressor, mask, lambda, [subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   glmodel - model number
    %   regressor - regressor name
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   lambda (optional) - regularization constant (defaults to 0, i.e. maximum likelihood estimate), or array of lambdas
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   dec - {nSubjects} struct array with decoded regressor for each subject:
    %      dec{s} - [nTRs x nVoxels x nLambdas] decoded regressor for each TR for each voxel in the mask, for subject s (optionally for each lambda)
    %
    % EXAMPLE:
    %   dec = ccnl_decode_regressor(exploration_expt(), 21, 'RU', 'rlpfc.nii', 1)
 

    if ~exist('lambda', 'var')
        lambda = 0;
    end

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

        % extract betas B and design matrix X
        cdir = pwd;
        cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
        if strcmp(mask_format, 'cor')
            B = spm_data_read(SPM.Vbeta, 'xyz', mask');
        else
            B = spm_data_read(SPM.Vbeta, mask);
        end
        cd(cdir);

        %X = SPM.xX.xKXs.X; %  use high-pass filtered & whitened design matrix
        X = SPM.xX.X; % use unfilterd, unwhitened design matrix

        % separate our regressor from the rest
        names = SPM.xX.name'; % regressor names
        which_reg = contains(names, regressor);

        assert(sum(which_reg) == length(EXPT.subject(subj).functional), 'Number of regressors that match is different from number of runs -- maybe add x and ^ prefix and suffix?');

        % separate X's and betas into matrices that do or don't have our regressor
        B_noreg = B(~which_reg, :);
        B_reg = B(which_reg, :);
        B_reg = repelem(B_reg, size(X, 1) / size(B_reg, 1), 1); % we need one for each TR b/c we're doing element-wise divison by b_RU
        X_noreg = X(:, ~which_reg);
        X_reg = X(:, which_reg);

        % extract activations 
        act = ccnl_get_activations(EXPT, glmodel, mask, subj);
        act = act{1};

        % decode regressor
        if length(lambda) == 1
            dec{s} = (act - X_noreg * B_noreg) .* B_reg ./ (B_reg.^2 + lambda);
        else
            numer = (act - X_noreg * B_noreg) .* B_reg;
            B_reg2 = B_reg.^2;
            for l = 1:length(lambda)
                dec{s}(:,:,l) = numer ./ (B_reg2 + lambda(l));
            end
        end
    end
