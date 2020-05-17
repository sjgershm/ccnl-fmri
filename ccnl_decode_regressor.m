function [dec, deconv, x, pmod] = ccnl_decode_regressor(EXPT, glmodel, regressor, mask, lambda, subjects, whiten_and_filter)

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
    %   [dec, deconv, x, pmod] = ccnl_decode_regressor(EXPT, glmodel, regressor, mask, [lambda,] [subjects,] [whiten_and_filter])
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
    %   whiten_and_filter (optional) -- whether to whiten & high-pass filter design matrix (default true)
    %
    % OUTPUTS:
    %   dec - {nSubjects} struct array with decoded regressor for each subject:
    %      dec{s} - [nTRs x nVoxels x nLambdas] decoded regressor for each TR for each voxel in the mask, for subject s (optionally for each lambda)
    %   deconv - {nSubjects} struct array with corresponding deconvolved regressor
    %      deconv{s} - [nOnsets x nVoxels x nLambdas] deconvolved decoded regressor for each event onset for each voxel in the mask, for subject s (optionally for each lambda)
    %   x - {nSubjects} struct array with original [nTRs x 1] regressor (convolved with HRF) for each subject; for comparison with dec
    %   pmod - {nSubjects} struct array with original [nOnsets x 1] parametric modulator for each subject; for comparison with deconv
    %
    % EXAMPLE:
    %   dec = ccnl_decode_regressor(exploration_expt(), 21, 'RU', 'rlpfc.nii', 1)
    %   
    %   [dec, deconv, x, pmod] = ccnl_decode_regressor(exploration_expt(), 45, 'masks/badre.nii', 1, 1)
    %   dec = mean(dec{1},2);
    %   deconv = mean(deconv{1}, 2);
    %   x = x{1};
    %   pmod = pmod{1};
    %   figure; plot([dec x]); legend({'decoded x HRF', 'original x HRF'}); xlabel('TR');
    %   figure; scatter(dec, x); xlabel('decoded x HRF'); ylabel('original x HRF');
    %   [r, p] = corr(dec, x)
    %   figure; plot([deconv pmod]); legend({'deconvolved', 'original'}); xlabel('TR');
    %   figure; scatter(deconv, pmod); pmodlabel('deconvolved'); ylabel('original');
    %   [r, p] = corr(deconv, pmod)
 

    if ~exist('lambda', 'var')
        lambda = 0;
    end

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('whiten_and_filter', 'var')
        whiten_and_filter = true;
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

        if whiten_and_filter
            X = SPM.xX.xKXs.X; %  use high-pass filtered & whitened design matrix
        else
            X = SPM.xX.X; % use unfilterd, unwhitened design matrix
        end

        % separate our regressor from the rest
        names = SPM.xX.name'; % regressor names
        which_reg = contains(names, regressor);
        disp(names(which_reg));

        assert(sum(which_reg) == length(EXPT.subject(subj).functional), 'Number of regressors that match is different from number of runs -- maybe add x and ^ prefix and suffix?');

        % separate X's and betas into matrices that do or don't have our regressor
        B_noreg = B(~which_reg, :);
        B_reg = B(which_reg, :);
        B_reg = repelem(B_reg, size(X, 1) / size(B_reg, 1), 1); % we need one for each TR b/c we're doing element-wise divison by b_RU
        X_noreg = X(:, ~which_reg);
        X_reg = X(:, which_reg);

        % extract activations 
        act = ccnl_get_activations(EXPT, glmodel, mask, subj, whiten_and_filter, whiten_and_filter);
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

        % also return regressor itself, for comparison
        x{s} = sum(X_reg, 2); % relies on random effects block diagonal design matrix
        assert(size(x{s}, 1) == size(dec{s},1), 'Decoded and original regressor have different lengths; contact Momchil...');

        %
        % deconvolve regressor:
        % fit a GLM with a beta series design matrix to the decoded regressor
        %

        fMRI_T     = SPM.xBF.T;
        fMRI_T0    = SPM.xBF.T0;

        % beta series design matrix -- from spm_fMRI_design
        Xx    = [];
        pmod{s} = [];
        for r = 1:length(SPM.nscan)
            %-Number of scans for this session
            %----------------------------------------------------------------------
            k = SPM.nscan(r);

            %-Get inputs, neuronal causes or stimulus functions U
            %----------------------------------------------------------------------
            U = spm_get_ons(SPM,r);

            % find regressor we want to decode
            found = false;
            for i = 1:length(U)
                j = find(contains(U(i).name, regressor));
                if ~isempty(j)
                    U = U(i);
                    assert(length(j) == 1, 'Found more than one matching regressor in the same run');
                    if j > 1
                        % it's a pmod
                        pmod{s} = [pmod{s}; U.P(j - 1).P];
                    else
                        % not a pmod => just add 1s
                        pmod{s} = [pmod{s}; ones(size(U.ons))];
                    end
                    found = true;
                    break;
                end
            end
            assert(found, 'Did not find regressor in U');

            % diaginalize it, so we can extract beta series
            % see spm_get_ons.m
            clear Udiag;
            for i = 1:length(U.ons)
                Udiag(i).name = {[U.name{1}, '_________', num2str(i)]};
                Udiag(i).ons = U.ons(i);
                Udiag(i).dur = U.dur(i);
                Udiag(i).P(1).name = 'none';
            end

            SPM.Sess(r).U = Udiag;
            Udiag = spm_get_ons(SPM,r);
            SPM.Sess(r).U = U;

            %-Convolve stimulus functions with basis functions
            %----------------------------------------------------------------------
            [X,Xn,Fc] = spm_Volterra(Udiag, SPM.xBF.bf, SPM.xBF.Volterra);

            %-Resample regressors at acquisition times (32 bin offset)
            %----------------------------------------------------------------------
            if ~isempty(X)
                X = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
            end

            %-Append into Xx
            %======================================================================
            Xx = blkdiag(Xx,X);
        end

        %-Design space and projector matrix [pseudoinverse] for WLS
        %---------------------------------------------------------
        % from spm_spm.m
        if whiten_and_filter
            xKXs        = spm_sp('Set', spm_filter(SPM.xX.K,SPM.xX.W * Xx));
        else
            warning('Using non-whitened and non-filtered design matrices violates the assumptions of the GLM');
            xKXs        = spm_sp('Set', Xx);
        end
        xKXs.X      = full(xKXs.X);
        pKX         = spm_sp('x-',xKXs);

        %-Weighted Least Squares estimation
        % =================================================
        for l = 1:length(lambda)
            beta        = pKX * dec{s}(:,:,l);
            deconv{s}(:,:,l) = beta;
        end

        assert(size(pmod{s}, 1) == size(deconv{s},1), 'Deconvolved regressor and pmod have different lengths; contact Momchil...');
    end
