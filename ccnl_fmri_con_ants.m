function ccnl_fmri_con_ants(EXPT,model,contrasts,subjects)
    
    % Construct contrasts and perform group-level analysis.
    
    % Note: this function should only be used after preprocessing data 
    % with ccnl_fmri_preproc_ants.m, and estimating subject-level GLMs
    % with ccnl_fmri_glm_ants.m. This function behaves exactly as
    % ccnl_fmri_con.m, but it normalizes the maps of beta coefficients and
    % t values using the transformations estimated during the preprocessing
    % phase.
    %
    % USAGE: ccnl_fmri_con_ants(EXPT,model,contrasts,[subjects])

    fprintf('Computing contrasts for GLM %d\n', model);
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end
    cdir = pwd;
    C = containers.Map({'+', '-'}, [1, -1]);
    spm('defaults','fmri');
    spm_jobman('initcfg');
    if isstr(contrasts); contrasts = {contrasts}; end
    
    [ccnl_fmri_dir, ~, ~] = fileparts(mfilename('fullpath')); 
    templatesdir_path = fullfile(ccnl_fmri_dir, 'templates');
    
    %% Construct contrasts
    for subj = subjects
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        convec = zeros(size(SPM.xX.name));
        for j = 1:length(contrasts)
            con = regexp(contrasts{j},'[-+]','split');
            sgn = regexp(contrasts{j},'[-+]','match');
            sgn = [{'+'}, sgn];
            contrasts{j}
            con
            sgn
            N = zeros(1,length(con));
            for c = 1:length(con)
                con{c} = strtrim(con{c});
                sgn{c} = strtrim(sgn{c});
                %found = false;
                ix = logical(zeros(1,length(SPM.xX.name)));
                for i = 1:length(SPM.xX.name)
                    if isequal(strfind(SPM.xX.name{i},[con{c},'*']), 1) || ~isempty(strfind(SPM.xX.name{i},[' ',con{c},'*'])) || ~isempty(strfind(SPM.xX.name{i},[' ',con{c},'^']))  || ~isempty(strfind(SPM.xX.name{i},['x',con{c},'^']))
                        convec(j,i) = C(sgn{c});
                        N(c) = N(c) + 1;
                        ix(i) = 1;
                    end
                end
                convec(j,ix) = convec(j,ix)/N(c);
                %found = true;
                %end
                %end
                %assert(found, ['Cannot find regressor ', con{c}]);
            end
            convec
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = contrasts{j};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = convec(j,:);
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
            
        end
        
        spm_jobman('run',matlabbatch);
        
        %% ANTs normalization of contrasts and T-map
        curdir_ants = pwd;
        
        disp('Normalizing contrasts and t-map using ANTs')
        S = EXPT.subject(subj);
        
       cd([EXPT.modeldir,'model',num2str(model),'/subj',num2str(subj)]);
        
        % normalize contrast betas
        con_files = dir('con*');
        for c = 1:length(con_files)
            command = ['antsApplyTransforms -d 3 -e 3 -v ' ...
                '-r ', templatesdir_path,'/TT_N27_SurfVol_2x2x2.nii ' ...
                '-i ', con_files(c).name, ' ' ...
                '-o w_', con_files(c).name, ' ' ...
                '-n Linear ' ...
                '-t ', S.datadir,'/w1Warp.nii.gz ' ...
                '-t ', S.datadir,'/w0GenericAffine.mat ' ...
                '-t [', S.datadir,'/anat2epi0GenericAffine.mat,1]'];
            system(command);
            command = ['mv w_', con_files(c).name,' ', con_files(c).name];
            system(command);
            
        end
        
        % normalize contrast t-map
        spmT_files = dir('spmT*');
        for c = 1:length(spmT_files)
            command = ['antsApplyTransforms -d 3 -e 3 -v ' ...
                '-r ', templatesdir_path,'/TT_N27_SurfVol_2x2x2.nii ' ...
                '-i ', spmT_files(c).name, ' ' ...
                '-o w_', spmT_files(c).name, ' ' ...
                '-n Linear ' ...
                '-t ', S.datadir,'/w1Warp.nii.gz ' ...
                '-t ', S.datadir,'/w0GenericAffine.mat ' ...
                '-t [', S.datadir,'/anat2epi0GenericAffine.mat,1]'];
            system(command);
            command = ['mv w_', spmT_files(c).name,' ', spmT_files(c).name];
            system(command);
        end
        
        cd(curdir_ants); % go back to whatever SPM was doing
        
    end
    
    cd(cdir);
    
    %% save contrasts
    modeldir0 = fullfile(EXPT.modeldir,['model',num2str(model)]);
    save(fullfile(modeldir0,'contrasts'),'contrasts','convec');
    
    %% Group-level analysis
    for j = 1:length(contrasts)
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(j)]);
        if isdir(modeldir); rmdir(modeldir,'s'); end
        mkdir(modeldir);
        cd(modeldir)
        job.dir{1} = modeldir;
        job.multi_cov = []; job.globalc.g_omit = 1;
        job.globalm.glonorm = 1; job.globalm.gmsca.gmsca_no = 1;
        job.masking.tm.tm_none = 1; job.masking.em{1} = [];
        
        % covariates - not working yet!!
        if isfield(EXPT,'cov')
            if ~isfield(EXPT.cov,'iCC')
                for i = 1:length(EXPT.cov); EXPT.cov(i).iCC = 1; end
            end
            if ~isfield(EXPT.cov,'iCFI')
                for i = 1:length(EXPT.cov); EXPT.cov(i).iCFI = 1; end
            end
            job.cov = EXPT.cov;
        else
            job.cov = [];
        end
        
        con = sprintf('con_%04d.nii',j);
        for s = 1:length(subjects)
            job.des.t1.scans{s} = fullfile(modeldir0,['subj',num2str(subjects(s))],con);
        end
        
        out = spm_run_factorial_design(job);
        load(out.spmmat{1});
        SPM = spm_spm(SPM);
        
        % write contrasts and t-maps
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
        matlabbatch{1}.spm.stats.con.delete = 1;
        for i = 1:length(SPM.xX.name)
            convec = zeros(size(SPM.xX.name)); convec(i) = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [contrasts{j},', ',SPM.xX.name{i}];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = convec';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        end
        
        spm_jobman('run',matlabbatch);
        
        cd(cdir);
    end
