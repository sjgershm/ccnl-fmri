function ccnl_mean(EXPT,subjects)
    
    % Write mean normalized structural.
    %
    % USAGE: ccnl_mean(EXPT,subjects)
    
    if nargin < 2; subjects = 1:length(EXPT.subject); end
    
    for subj = 1:length(subjects)
        i = subjects(subj);
        P{subj} = fullfile(EXPT.subject(i).datadir,'wBrain.nii');
    end
    
    spm_mean(P);
    movefile('mean.nii',fullfile(EXPT.modeldir,'mean.nii'));