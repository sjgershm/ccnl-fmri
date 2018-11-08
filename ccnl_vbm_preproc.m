function ccnl_vbm_preproc(EXPT,subjects,def)
   
    % Preprocess structural fMRI images for voxel-based morphometry (VBM).
    % Does tissue segmentation, spatial normalization (warping), and modulation.
    %
    % USAGE: ccnl_vbm_preproc(EXPT,[subjects],[def])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   subjects (optional) - vector of subject numbers (default: all subjects)
    %   def (optional) - parameter defaults (def = spm_get_defaults)
    %
    % OUTPUTS:
    %   Creates tissue maps such as mwc1struct.nii
    %
    % Momchil Tomov, Nov 2018
    
    % Setup
    if nargin < 3; def = spm_get_defaults; end % SPM defaults
    if nargin < 2; subjects = 1:length(EXPT.subject); end
    
    % Loop over subjects
    for subj = subjects
        
        S = EXPT.subject(subj); % subject structure
        for i = 1:length(S.functional); S.functional{i} = fullfile(S.datadir,S.functional{i}); end
        S.structural = fullfile(S.datadir,S.structural);
        
        %% Segment structural image
        job = [];
        ngaus  = [1 1 2 3 4 2];
        warped = [1 1 1 0 0 0];
        modulated = [1 1 1 0 0 0];
        for c = 1:6 % tissue classes
            job.tissue(c).tpm = {fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
            job.tissue(c).ngaus = ngaus(c);
            job.tissue(c).native = [0 0];
            job.tissue(c).warped = [warped(c) modulated(c)];
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
    end
