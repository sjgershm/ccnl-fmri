function ccnl_mean(EXPT,subjects)
    
    % Write mean normalized structural.
    %
    % USAGE: ccnl_mean(EXPT,subjects)
    
    if nargin < 2; subjects = 1:length(EXPT.subject); end
    
    for subj = 1:length(subjects)
        i = subjects(subj);
        pth = fileparts(EXPT.subject(i).structural);
        P{subj} = fullfile(pth,'wBrain.nii');
    end
    
    spm_mean(P);