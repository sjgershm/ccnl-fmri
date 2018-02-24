function ccnl_fmri_glm(EXPT,model,subjects)
    
    % Estimate subject-level GLM.
    %
    % USAGE: ccnl_fmri_glm(EXPT,model,[subjects])
    
    cdir = pwd;
    if nargin < 3; subjects = 1:length(EXPT.subject); end
    
    % create models folder if none exists
    if ~isdir(EXPT.modeldir); mkdir(EXPT.modeldir); end
    if ~isdir(fullfile(EXPT.modeldir,['model',num2str(model)])); mkdir(fullfile(EXPT.modeldir,['model',num2str(model)])); end
    
    % generic design specification
    def = spm_get_defaults;
    job.timing.RT = EXPT.TR;
    job.timing.units = 'secs';
    job.timing.fmri_t = def.stats.fmri.t;
    job.timing.fmri_t0 = def.stats.fmri.t0;
    job.volt = 1;
    job.fact = [];
    job.cvi = def.stats.fmri.cvi;
    job.global = 'None';
    job.mthresh = 0.4;
    
    % hrf specification; use canonical hrf by default
    if isfield(EXPT,'bases')
        job.bases = bases;
    else
        job.bases.hrf.derivs = [0 0];
    end
    
    job0 = job;
    
    % subject-specific design specification
    for subj = subjects
        
        job = job0;
        
        % overwrite existing files
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        if isdir(modeldir); rmdir(modeldir,'s'); end
        mkdir(modeldir);
        
        job.dir{1} = modeldir;
        S = EXPT.subject(subj);
        job.mask{1} = fullfile(S.datadir,'wBrain.nii');
        for i = 1:length(S.functional)
            multi = EXPT.create_multi(model,subj,i);
            save(fullfile(modeldir,['multi',num2str(i)]),'-struct','multi');
            job.sess(i).hpf = def.stats.fmri.hpf;   % high-pass filter
            job.sess(i).scans{1} = fullfile(S.datadir,['swu',S.functional{i}]);
            job.sess(i).multi{1} = fullfile(modeldir,['multi',num2str(i)]);
            job.sess(i).cond = struct('name',{},'onset',{},'duration',{},'tmod',{},'pmod',{},'orth',{});
            job.sess(i).regress = [];
            job.sess(i).multi_reg{1} = spm_file(fullfile(S.datadir,S.functional{i}),'prefix','rp_','ext','txt');    % motion regressors from realignment
        end
        
        cd(modeldir);

        job = spm_run_fmri_spec(job);
        load(job.spmmat{1});
        spm_spm(SPM);

        cd(cdir);
    end
