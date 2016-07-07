function ccnl_fmri_preproc(EXPT,subjects)
    
    % Preprocess fMRI data. This function does the following: (1)
    % realignment and unwarping, (2) segmentation, (3) coregistration of functionals to
    % structural, (4) normalization, (5) smoothing.
    %
    % USAGE: fmri_preproc(EXPT,subjects)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   subjects (optional) - vector of subject numbers (default: all subjects)
    %
    % OUTPUTS:
    %   Realignment and coregistration induce changes in the headers of the nifti files.
    %   Unwarping writes out new files with 'u' prefix.
    %   Normalization writes out new files with 'w' prefix.
    %   Smoothing writes out new files with 's' prefix.
    %
    % Sam Gershman, June 2016
    
    % Setup
    def = spm_get_defaults; % SPM defaults
    if nargin < 2; subjects = 1:length(EXPT.subject); end
    
    % Loop over subjects
    for subj = subjects
        
        S = EXPT.subject(subj); % subject structure
        for i = 1:length(S.functional); S.functional{i} = fullfile(S.datadir,S.functional{i}); end
        S.structural = fullfile(S.datadir,S.structural);
        
        %% Realign functional images; write movement parameters to rp*.txt
        spm_realign(S.functional);
        
        %% Unwarp functional images
        for i = 1:length(S.functional)
            ds = spm_uw_estimate(S.functional{i});  % estimate deformation field
            spm_uw_apply(ds);                       % reslice unwarped functionals (u* prefix)
        end
        
        %% Segment structural image
        job = [];
        ngaus  = [1 1 2 3 4 2];
        native = [1 1 1 0 0 0];
        for c = 1:6 % tissue classes
            job.tissue(c).tpm = {fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
            job.tissue(c).ngaus = ngaus(c);
            job.tissue(c).native = [native(c) 0];
            job.tissue(c).warped = [0 0];
        end
        job.channel.biasreg = 0.001;            % light regularization
        job.channel.biasfwhm = 60;              % FWHM of Gaussian smoothness of bias
        job.channel.write = [0 1];              % save bias-corrected image
        job.warp.mrf = 1;                       % MRF cleanup parameter
        job.warp.cleanup = 1;                   % light MRF cleanup
        job.warp.fwhm = 0;
        job.warp.affreg = 'mni';
        job.warp.reg = [0 0.001 0.5 0.05 0.2];
        job.warp.samp = 3;
        job.warp.write = [0 1];
        job.channel.vols{1} = S.structural;
        vout = spm_preproc_run(job);
        
        %% Write skull-stripped bias-corrected structural (Brain.nii)
        for i = 1:3
            Vi(i) = spm_vol(vout.tiss(i).c{1});
        end
        Vi(4) = spm_vol(vout.channel.biascorr{1});
        f = '(i1 + i2 + i3) .* i4';
        Vstruc = spm_imcalc(Vi,fullfile(S.datadir,'Brain.nii'),f);
        
        %% Coregister structural and mean functional
        disp('Coregistering structural and mean functional...')
        job = [];
        job.eoptions = def.coreg.estimate;
        job.ref{1} = Vstruc.fname;
        job.source{1} = spm_file(S.functional{1},'prefix','meanu','number','');
        job.other = spm_file(S.functional,'prefix','u');
        out = spm_run_coreg(job);
        
        %% Normalize functionals using forward deformation parameters from segmentation
        disp('Normalizing functionals...')
        job = [];
        job.woptions = def.normalise.write;
        job.subj.resample = out.cfiles;
        job.subj.def{1} = vout.fordef{1};
        out = spm_run_norm(job);
        
        %% Smooth normalized functionals
        disp('Smoothing functionals...')
        job = def.smooth;
        job.data = out.files;
        job.dtype = 0;
        job.im = 0;             % no implicit masking
        spm_run_smooth(job);
        
        %% Normalize skull-stripped bias-corrected structural
        disp('Normalizing structural...')
        job = [];
        job.woptions = def.normalise.write;
        job.woptions.vox = [1 1 1];
        job.subj.resample{1} = Vstruc.fname;
        job.subj.def{1} = vout.fordef{1};
        spm_run_norm(job);
        
        %% Write mean functional image (only used for checking coregistration)
        spm_mean(spm_file(S.functional{1},'prefix','wu'));
        movefile('mean.nii',fullfile(S.datadir,'mean.nii'));
        
    end