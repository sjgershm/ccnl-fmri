function ccnl_view(model,contrast)
    
    % View group-level T-map. If 'mean.nii' (average normalized structural) exists, will use that as overlay.
    %
    % USAGE: ccnl_view(model,contrast)
    %
    % Requires bspmview (http://www.bobspunt.com/bspmview/).
    
    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
    load(fullfile(modeldir,'contrasts'));
    ix = find(strcmp(contrasts,contrast));
    if isempty(ix)
        error('Contrast not found');
    end
    spmT = fullfile(EXPT.modeldir,['model',num2str(model)],sprintf('spmT_%04d.nii',ix));
    struc = fullfile(EXPT.modeldir,'mean.nii');
    if ~exist(struc,'file')
        bspmview(spmT);
    else
        bspmview(spmT,struc);
    end