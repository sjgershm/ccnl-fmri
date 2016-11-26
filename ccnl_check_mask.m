function ccnl_check_mask(EXPT,model,contrast)
    
    % Check co-registration of group-level structural and mask images.
    % If 'mean.nii' (average normalized structural) doesn't exist,
    % one will be created.
    %
    % USAGE: ccnl_check_mask(EXPT,model,contrast)

    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
    load(fullfile(modeldir,'contrasts'));
    ix = find(strcmp(contrasts,contrast));
    if isempty(ix)
        error('Contrast not found');
    end
    mask = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(ix)],'mask.nii');
    struc = fullfile(EXPT.modeldir,'mean.nii');
    if ~exist(struc,'file')
        ccnl_mean(EXPT);
    end

    P{1} = mask;
    P{2} = struc;
    spm_check_registration(char(P));
