function ccnl_view(EXPT,model,contrast)
    
    % View group-level T-map. If 'mean.nii' (average normalized structural)
    % exists, one will be created for use as an overlay.
    %
    % USAGE: ccnl_view(EXPT,model,contrast)
    %
    % Requires bspmview (http://www.bobspunt.com/bspmview/).
    
    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
    load(fullfile(modeldir,'contrasts'));
    ix = find(strcmp(contrasts,contrast));
    if isempty(ix)
        error('Contrast not found');
    end
    spmT = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(ix)],'spmT_0001.nii');
    struc = fullfile(EXPT.modeldir,'mean.nii');
    if ~exist(struc,'file')
        ccnl_mean(EXPT);
    end
    
    bspmview(spmT,struc);