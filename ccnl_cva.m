function CVA = ccnl_cva(EXPT,model,cons,mask,subjects)

    % Canonical variate analysis (CVA). Test how many mutually independent relations there are
    % between given contrasts and multivariate pattern of activity in masked region.
    %
    % See spm_cva.m and spm_cva_ui.m for more info.
    % Also Darlington et al. 1973 for an nice explanation of CVA.
    %
    % USAGE:
    %   CVA = ccnl_cva(EXPT, model, cons, mask [, subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number (matters only if whitening, otherwise any one will do)
    %   cons - cell array of contrast names to use as features
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   CVA{s} - CVA result struct (see spm_cva.m), most importantly:
    %     .p  - [nCons] p-value array for testing that the number of relations is >= i (i.e. for rejecting the null that it's < i)
    % 
    % EXAMPLE:
    %   CVA = ccnl_cva(vgdl_expt(), 1, {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda'}, 'sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii', 1)

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    % find contrast indices
    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
    load(fullfile(modeldir,'contrasts'));
    ix = find(ismember(contrasts,cons));

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isempty(Vmask) || isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');

        % get SPM-friendly MNI voxel coordinates
        if strcmp(mask_format, 'cor')
            XYZ = mask;
        else
            assert(ismember(mask_format, {'mask', 'inds_or_binary'}));
            [x y z] = ind2sub(SPM.xY.VY(1).dim, mask);
            XYZ = [x y z];
        end
        XYZ = cor2mni(XYZ, SPM.xVol.M); 
        XYZ = XYZ';

        %-Contrast specification
        %------------------------------------------------------------------
        con    = {SPM.xCon(ix).name}; % momchil note: there is a bug in SPM; it doesn't handle multiple contrasts properly
        c      = [SPM.xCon(ix).c];
        c      = full(c);
        
        % extract voxel activations
        Y = ccnl_get_activations(EXPT, model, mask, subj, true, false); % note they don't high-pass filter in spm_cva_ui.m
        Y = Y{1};

        %-Remove serial correlations and get design (note X := W*X)
        %------------------------------------------------------------------
        X   = SPM.xX.xKXs.X;
        
        %-Null-space
        %------------------------------------------------------------------
        X0  = [];
        try, X0 = [X0 blkdiag(SPM.xX.K.X0)]; end          %-drift terms
        try, X0 = [X0 spm_detrend(SPM.xGX.gSF)]; end      %-global estimate
        
        %-Canonical Variate Analysis
        %==================================================================
        U     = spm_mvb_U(Y,'compact',spm_svd([X0, X-X*c*pinv(c)]),XYZ);
        CVA{s}   = spm_cva(Y, X, X0, c, U);


        %-Assemble results
        %------------------------------------------------------------------
        CVA{s}.contrast = con;                      %-contrast name
        CVA{s}.XYZ      = XYZ;                      %-locations of voxels (mm)
        CVA{s}.U        = U;                        %-dimension reduction (SVD)
        CVA{s}.c        = c;
    end   
