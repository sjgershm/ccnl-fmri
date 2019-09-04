function ccnl_fmri_preproc_ants(EXPT,subjects,def)

% Preprocess fMRI data using a combination of SPM and the Advanced
% Normalization Tools (ANTs). The latter is a toolbox that achieves
% superior registration of participants brain to a template, using
% diffeomorphic registration based on prior knowledge maps.
%
% Note that this function requires ANTs to be available to the OS. 
%
% This function does the following: (1) realignment and unwarping 
% using SPM, (2) smoothing of the functionals using SPM, (3) normalization 
% of structural to the TT_N27_SurfVol template using ANTs, and 
% (4) estimation of transformation from the structural to anatomical
% using ANTs.
%
% If you use this pipeline, you will be set up to do your analyses in 
% participant's native space (i.e., on their non-normalized functionals). 
% You will have to do this using ccnl_fmri_glm_ants.m. Then you will use 
% ccnl_fmri_con.m to normalize the statistical maps before you can perform 
% inferential statistics.
%
% USAGE: fmri_preproc_ants(EXPT,[subjects],[def])
%
% INPUTS:
%   EXPT - experiment structure
%   subjects (optional) - vector of subject numbers (default: all subjects)
%   def (optional) - parameter defaults (def = spm_get_defaults)
%
% OUTPUTS:
%   Realignment induces changes in the headers of the nifti files.
%   Unwarping writes out new files with 'u' prefix.
%   Normalization writes out new files with 'w' prefix.
%   Smoothing writes out new files with 's' prefix.
%
% Wouter Kool, July 2019 (based on code by Sam Gershman, June 2016)

% Setup
if nargin < 3; def = spm_get_defaults; end % SPM defaults
if nargin < 2; subjects = 1:length(EXPT.subject); end

[ccnl_fmri_dir, ~, ~] = fileparts(mfilename('fullpath')); 
templatesdir_path = fullfile(ccnl_fmri_dir, 'templates');

% Loop over subjects
for subj = subjects
    
    S = EXPT.subject(subj); % subject structure
    for i = 1:length(S.functional); S.functional{i} = fullfile(S.datadir,S.functional{i}); end
    S.structural = fullfile(S.datadir,S.structural);
    
    %% Realign functional images; write movement parameters to rp*.txt
    spm_realign(S.functional);    % estimate alignment
 
    %% Unwarp functional images
    for i = 1:length(S.functional)
        ds = spm_uw_estimate(S.functional{i});                   % estimate deformation field
        spm_uw_apply(ds);                                     % write unwarped functionals (u* prefix)
    end
    
    %% Reslice functional images
    spm_reslice(spm_file(S.functional,'prefix','u'));      % write realigned functionals (r* prefix)
    
    %% Smooth  functionals
    disp('Smoothing functionals...')
        
    job = def.smooth;
    job.data = spm_select('expand',spm_file(S.functional,'prefix','ru'));
    job.dtype = 0;
    job.im = 0;             % no implicit masking
    spm_run_smooth(job);
    
    %% Skull strip structural image using ANTs
    currentdir = pwd;
    cd(S.datadir); % getting ready for ANTs
    
    disp('Skull-stripping structural...')
    
    command = ['antsBrainExtraction.sh -d 3 -a struct.nii -o brain ' ...
        '-e ',templatesdir_path,'T_template0.nii.gz ' ...
        '-m ',templatesdir_path,'T_template0_BrainCerebellumProbabilityMask.nii.gz ' ...
        '-f ',templatesdir_path,'T_template0_BrainCerebellumRegistrationMask.nii.gz'];
    system(command);
    system('mv brainBrainExtractionBrain.nii.gz Brain.nii.gz');
    
    %% Normalize skull-stripped brain to template using ANTs
    disp('Normalizing skull-stripped structural..')
    
    command = ['antsRegistrationSyN.sh -d 3 -m Brain.nii.gz ' ...
        '-f ',templatesdir_path,'TT_N27_SurfVol.nii -o w -n 8'];
    system(command);
    system('mv wWarped.nii.gz wBrain.nii.gz');
    system('gunzip -f wBrain.nii.gz');
    
    %% Coregister structural to mean functional of first run using ANTs
    disp('Coregister structural to mean functional...')
    
    command = ['antsRegistrationSyN.sh -f ',spm_file(S.functional{1},'prefix','meanu','number',''),' -m Brain.nii.gz -o anat2epi -t r'];
    system(command);
    
    system('mv anat2epiWarped.nii.gz mask.nii.gz');
    system('gunzip -f mask.nii.gz');
    
    cd(currentdir);
    
end
