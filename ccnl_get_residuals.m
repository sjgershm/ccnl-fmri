function residuals = ccnl_get_residuals(EXPT,model,mask,subjects,whiten,filter)
    
    % Extract residuals coefficients from a mask.
    % Caution: don't use for too many voxels
    %
    % USAGE: residuals = ccnl_get_residuals(EXPT,model,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %   whiten (optional) - whether to whiten data like SPM does (default false)
    %   filter (optional) - whether to high-pass filter data like SPM does (default false)
    %
    % OUTPUTS:
    %   residuals{s} - [nScans x nVoxels] residuals for subject s
    %
    % Momchil Tomov, Aug 2018
    
    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('whiten', 'var')
        whiten = false;
    end

    if ~exist('filter', 'var')
        filter = false;
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    for s = 1:length(subjects)
        subj = subjects(s);

        Y = ccnl_get_activations(EXPT, model, mask, subj, whiten, filter);
        Y = Y{1};

        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isempty(Vmask) || isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');

        X = SPM.xX.X;
        if whiten
            X = SPM.xX.W*X;
        end
        if filter
            X = spm_filter(SPM.xX.K,X);
        end
        Xs = spm_sp('Set', X);
        Xs.X = full(Xs.X); % see spm_spm.m

        % extract residuals
        res = spm_sp('r',Xs,Y);

        % sanity checks
        %
        if ~whiten & ~filter
            assert(immse(SPM.xX.X, Xs.X) < 1e-10);
        end
        if whiten & filter
            assert(immse(SPM.xX.xKXs.X, Xs.X) < 1e-10);

            % compare to sum of squared residuals computed by SPM
            ResSS = sum(res.^2); %-Residual SSQ
            V = spm_vol(fullfile(modeldir,'ResMS.nii'));
            if strcmp(mask_format, 'cor')
                ResMS = spm_data_read(V, 'xyz', mask');
            else
                ResMS = spm_data_read(V, mask);
            end
            assert(immse(ResSS / SPM.xX.trRV, ResMS) < 1e-9, ['Computed residuals don''t match ResMS.nii for subject', num2str(subj), ' -- maybe mask has out-of-brain voxels?']);  % ResMS = ResSS scaled by tr(RV)
        end

        residuals{s} = res;

        fprintf('Computed residuals for subject %d\n', subj);
    end
