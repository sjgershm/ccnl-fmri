function ccnl_fmri_con(EXPT,model,contrasts,subjects)
    
    % Construct contrasts and perform group-level analysis.
    %
    % USAGE: ccnl_fmri_con(EXPT,model,contrasts,[subjects])
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end
    cdir = pwd;
    C = [1 -1];
    spm('defaults','fmri');
    spm_jobman('initcfg');
    
    %% Construct contrasts
    for subj = subjects
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        for j = 1:length(contrasts)
            con = regexp(contrasts{j},'-','split');
            convec = zeros(size(SPM.xX.name));
            for c = 1:length(con)
                con{c} = strtrim(con{c});
                for i = 1:length(SPM.xX.name)
                    if ~isempty(strfind(SPM.xX.name{i},[con{c},'*'])) || ~isempty(strfind(SPM.xX.name{i},[con{c},'^']))
                        convec(i) = C(c);
                    end
                end
            end
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = contrasts{j};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = convec;
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        end
        
        spm_jobman('run',matlabbatch);
    end
    
    cd(cdir);
    
    %% Group-level analysis
    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
    save(fullfile(modeldir,'contrasts'),'contrasts');
    cd(modeldir)
    job.dir{1} = modeldir;
    job.multi_cov = []; job.cov = []; job.globalc.g_omit = 1;
    job.globalm.glonorm = 1; job.globalm.gmsca.gmsca_no = 1;
    job.masking.tm.tm_none = 1; job.masking.em{1} = [];
    for j = 1:length(contrasts)
        con = sprintf('con_%04d.nii',j);
        for s = 1:length(subjects)
            job.des.t1.scans{s} = fullfile(modeldir,['subj',num2str(subjects(s))],con);
        end
    end
    out = spm_run_factorial_design(job);
    load(out.spmmat{1});
    SPM = spm_spm(SPM);
    
    % write contrasts and t-maps
    matlabbatch = [];
    matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
    matlabbatch{1}.spm.stats.con.delete = 1;
    for j = 1:length(contrasts)
        convec = zeros(size(SPM.xX.name));
        ix = strcmp(SPM.xX.name,'mean');
        convec(ix) = 1;
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = contrasts{j};
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = convec;
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
    end
    
    spm_jobman('run',matlabbatch);
    
    cd(cdir);