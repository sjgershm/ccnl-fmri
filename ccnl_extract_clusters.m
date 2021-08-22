function [V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, model, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent, df)
    %
    % Given a contrast, extract all the activation clusters from the t-map 
    % after cluster FWE correction. Code is copy-pasted and adapted from
    % bspmview version 20161108. 
    % As a sanity check prints the results table -- should be the same as
    % the one from bspmview.
    %
    % USAGE:
    %   [V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, model, contrast, p, direct, alpha, Dis, Num)
    %
    % EXAMPLES:
    %   [V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(exploration_expt(), 1, 'RR - SS')
    %   [V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(exploration_expt(), 1, 'RR - SS', 0.001, '+/-', 0.05, 20, 3)
    %
    % INPUT:
    %   EXPT = experiment structure, e.g. exploration_expt()
    %   model = GLM number, e.g. 1
    %   contrast = contrast, e.g. 'RR - SS'
    %   p (optional) = p-value threshold; defaults to 0.001
    %   direct (optional) = sign of the activations; should be one of '+', '-', or '+/-'; defaults to '+/-'
    %   alpha (optional) = significance level for cluster FWE correction, default 0.05 in bspmview
    %   Dis (optional) = separation for cluster maxima, default 20 in bspmview
    %   Num (optional) = numpeaks for cluster maxima, default 1 here (3 in bspmview)
    %   clusterFWEcorrect (optional) = whether to perform cluster FWE correction (defaults to true)
    %   extent (optional) = minimum cluster extent (overrides clusterFWEcorrect; not used by default)
    %   df (optional) = degrees of freedom in t-test (# of subjects - 1)
    %
    % OUTPUT:
    %   V = SPM volume of the t-map, with the filename changed so we don't overwrite it accidentally
    %   Y = the actual t-map
    %   C = volume with cluster size for each voxel
    %   CI = volume with cluster index for each voxel (i.e. the actual clusters) 
    %   region = labels for the peak voxels based on AAL2 atlas
    %   extent = size of each cluster
    %   stat = t-statistic for the peak voxels
    %   mni = MNI coordinates of peak voxels
    %   cor = coordinates of peak voxel in native space (can be used as indices in C and CI)
    %   results_table = what Save Results Table in bspmview would spit out 
    %   spmT = path to the t-map .nii file
    %
    % Momchil Tomov, Sep 2018
    % 

    if ~exist('p', 'var')
        p = 0.001;
    end
    if ~exist('direct', 'var')
        direct = '+/-';
    end
    if ~exist('alpha', 'var')
        alpha = 0.05;
    end
    if ~exist('Dis', 'var')
        Dis = 20;
    end
    if ~exist('Num', 'var')
        Num = 1; % default 3 in bspmview
    end
    if ~exist('clusterFWEcorrect', 'var')
        clusterFWEcorrect = true;
    end
    if ~exist('extent', 'var')
        extent = [];
    end

    % only used for sanity check
    atlas_name = 'AAL2';
    [curdir, ~, ~] = fileparts(mfilename('fullpath')); 
    atlas_dirpath = fullfile(curdir, 'atlases');

    assert(ismember(direct, {'+/-', '+', '-'}));

    if ischar(EXPT)
        spmT = EXPT; % HACK TODO FIXME -- this is so I can pass a tmap directly if I want to
    else
        % find the contrast
        %
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)]);
        load(fullfile(modeldir,'contrasts'));
        ix = find(strcmp(contrasts,contrast));
        if isempty(ix)
            error('Contrast not found');
        end
        spmT = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(ix)],'spmT_0001.nii');
    end

    if ~exist('df', 'var')
        SPM_path = fullfile(EXPT.modeldir,['model',num2str(model)],['con',num2str(ix)],'SPM.mat');
        load(SPM_path);
        % get degrees of freedom from SPM file
        df = SPM.nscan - 1; 
    end

    fprintf('bspm_extract_clusters(%s, %f, %s, %f, %d, %d, %d, %s, %s, %d, %d)\n', spmT, p, direct, alpha, Dis, Num, df, atlas_dirpath, atlas_name, clusterFWEcorrect, extent);


    % extract the clusters
    %
    [C, CI, region, extent, stat, mni, cor, results_table] = bspm_extract_clusters(spmT, p, direct, alpha, Dis, Num, df, atlas_dirpath, atlas_name, clusterFWEcorrect, extent);

    V = spm_vol(spmT);
    Y = spm_read_vols(V);
    V.fname = []; % <-- change immediately so we don't overwrite it accidentally

    disp(results_table);
end





function [C, CI, region, extent, stat, mni, cor, results_table] = bspm_extract_clusters(tmap_filename, p, direct, alpha, Dis, Num, df, atlas_dirpath, atlas_name, clusterFWEcorrect, extent)
    %
    % Ripped and refactored from bspmview.
    % Extract clusters and peak voxels from a t-map contrast after cluster FWE
    % correction. Exactly the same as bspmview -- uses the same functions. As a
    % bonus, even spits out a results table.
    %
    % INPUT:
    %   tmap_filename = path to .nii file with the contrast t-map, e.g.
    %                   '../neural/model154/con10/spmT_0001.nii'
    %   p = p-value threshold for the individual voxels, e.g. 0.001
    %   direct = sign of the activations to look at; should be one of +, -, or
    %            +/-
    %   alpha = significance level for cluster FWE correction, default 0.05 in bspmview
    %   Dis = separation for cluster maxima, default 20 in bspmview
    %   Num = numpeaks for cluster maxima, default 3 in bspmview
    %   df = degrees of freedom for t-tests = # subjects - 1
    %   atlas_dirpath = path to bspmview atlas directory
    %   atlas_name = bspmview atlas name
    %   clusterFWEcorrect = whether to do cluster FWE correction
    %   extent = min cluster extent (overrides cluster FWE corr) 
    % 
    % OUTPUT:
    %   C = volume with cluster size for each voxel
    %   CI = volume with cluster index for each voxel <-- that's the name of the
    %        game; 
    %   region = labels for the peak voxels
    %   extent = size of the cluster
    %   stat = t-statistic for the peak voxels
    %   mni = MNI coordinates of peak voxels
    %   cor = coordinates of peak voxel in native space (can be used as indices
    %         in the C and CI volumes)
    %   results_table = what Save Results Table in bspmview would spit out 
    % 

    %tmap_filename = '../neural/model154/con10/spmT_0001.nii';
    %p = 0.001; % p-value threshold for individual voxels
    %direct = '+/-'; % positive or negative activations
    %alpha = 0.05;
    %Dis = 20;
    %Num = 1;


    assert(ismember(direct, {'+', '-', '+/-'}));


    % get cluster extent threshold
    %
    if isempty(extent)
        extent_thresh = bspm_cluster_correct(tmap_filename, df, direct, p, alpha); % FWE; also supports FDR as second argument
        if isinf(extent_thresh) || ~clusterFWEcorrect 
            warning('No voxels found after cluster FWE correction. Setting extent threshold = 5.');
            extent_thresh = 5;
        end
    else
        extent_thresh = extent;
    end
    fprintf('extent threshold = %d\n', extent_thresh);


    % get cluster indices
    %
    thresh = spm_invTcdf(1-p, df);
    
    p
    df
    thresh
    
    V = spm_vol(tmap_filename);
    Y = spm_read_vols(V);
    Y(isnan(Y)) = 0;
    [clustsize, clustidx]   = bspm_getclustidx(Y, thresh, extent_thresh);
    di = strcmpi({'+' '-' '+/-'}, direct);
    C = clustsize(di, :);
    CI= clustidx(di, :);


    % set some input params
    %
    M           = V.mat;         %-voxels to mm matrix
    DIM         = V.dim';
    VOX         = abs(diag(M(:,1:3)));
    [x,y,z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ0         = [x(:)';y(:)';z(:)'];
    RCP         = XYZ0;
    RCP(4,:)    = 1;
    XYZmm0       = M(1:3,:)*RCP;


    % set thresh
    %
    idx = find(C > 0);
    if find(di) == 2
        Z     = abs(Y(idx));
    else
        Z     = Y(idx);
    end
    Nunique   = length(unique(Z));
    XYZ       = XYZ0(:,idx);
    XYZmm     = XYZmm0(:,idx);


    % set maxima
    %
    switch char(direct)
        case {'+', '-'}
            LOCMAX      = bspm_getmaxima(Z, XYZ, M, Dis, Num);
            if strcmp(direct, '-')
                LOCMAX(:,2) = -LOCMAX(:,2);
            end
        otherwise
            POS         = bspm_getmaxima(Z, XYZ, M, Dis, Num);
            NEG         = bspm_getmaxima(Z*-1, XYZ, M, Dis, Num);
            if ~isempty(NEG), NEG(:,2) = NEG(:,2)*-1; end
            LOCMAX      = [POS; NEG];
    end


    % get labels
    %
    [atlaslabels, atlas] = bspm_setatlas(tmap_filename, atlas_dirpath, atlas_name);
    LABELS = bspm_getregionnames(LOCMAX(:,3:5)', atlaslabels, atlas, XYZmm0);

    % generate resultstable
    %
    results_table = [cell(size(LABELS)) LABELS num2cell(LOCMAX)];
    results_table{1,1} = 'Positive';
    if any(LOCMAX(:,2)<0)
        tmpidx = find(LOCMAX(:,2)<0);
        results_table{tmpidx(1),1} = 'Negative';
    end

    % set outputs
    %
    C = reshape(C, size(Y));
    CI = reshape(CI, size(Y));
    region = LABELS;
    extent = LOCMAX(:, 1);
    stat = LOCMAX(:, 2);
    mni = LOCMAX(:, 3:5);
    cor = mni2cor(mni, V.mat);

    %% for sanity checks w/ the real bspmview
    %
    %{

    LOCMAX_orig = LOCMAX;

    load('~/Downloads/bspm_setmaxima.mat');

    assert(immse(LOCMAX_orig, LOCMAX) < 1e-5);


    C_orig = C;
    CI_orig = CI;


    load('~/Downloads/bspm_correct.mat');
    disp(unique(C));
    assert(isequal(C_orig, C));

    [st.ol.C0, st.ol.C0IDX] = bspm_getclustidx(st.ol.Y, thresh, extent);
    C = st.ol.C0(di,:);
    CI= st.ol.C0IDX(di, :);
    disp(unique(C));
    disp(unique(CI));
    assert(isequal(C_orig, C));
    assert(isequal(CI_orig, CI));


    load('~/Downloads/bspm_od.mat')
    [clustsize, clustidx]   = bspm_getclustidx(od, thresh, extent);
    C = clustsize(di, :);
    CI= clustidx(di, :);
    disp(unique(C));
    disp(unique(CI));
    assert(isequal(C_orig, C));
    assert(isequal(CI_orig, CI));

    % '../neural/model154/con10/spmT_0001.nii'
    %}
end





function [fwek, fdrk] = bspm_cluster_correct(im,df,direct,u,alpha)
    % Cluster correction function copy-pasted from bspmview, with some modifications to 
    % make it work standalone.
    %
    % INPUT:
    % im = Path to .nii file with t-map, e.g. '../neural/model154/con10/spmT_0001.nii'
    % df = degrees of freedom = # of subjects - 1 (for group-level contrast)
    % direct = direction of t-values, should be one of {'+' '-' '+/-'}
    % u = p-value threshold (exclude all voxels with greater p-values), default 0.001
    % alpha = significance level of the cluster correction (default 0.05)
    %
    % OUTPUT:
    % fwek = Cluster FWE extent threshold
    % fdrk = Cluster FDR extent threshold
    % 
    % THIS IS A MODIFICATION OF A FUNCTION BY BSPMVIEW. THIS IS THE ORIGINAL DOCUMENTATION
    %
    % CLUSTER_CORRECT Computer extent for cluster-level correction
    % USAGE: [k info] = cluster_correct(im,u,alpha,range)
    %
    %
    % THIS IS A MODIFICATION OF A FUNCTION BY DRS. THOMAS NICHOLS AND MARKO
    % WILKE, CorrClusTh.m. ORIGINAL DOCUMENTATION PASTED BELOW:
    %
    % Find the corrected cluster size threshold for a given alpha
    % function [k,Pc] =CorrClusTh(SPM,u,alpha,guess)
    % SPM   - SPM data structure
    % u     - Cluster defining threshold
    %         If less than zero, u is taken to be uncorrected P-value
    % alpha - FWE-corrected level (defaults to 0.05)
    % guess - Set to NaN to use a Newton-Rhapson search (default)
    %         Or provide a explicit list (e.g. 1:1000) of cluster sizes to
    %         search over.
    %         If guess is a (non-NaN) scalar nothing happens, except the the
    %         corrected P-value of guess is printed.
    %
    % Finds the corrected cluster size (spatial extent) threshold for a given
    % cluster defining threshold u and FWE-corrected level alpha.
    %
    %_________________________________________________________________________
    % $Id: CorrClusTh.m,v 1.12 2008/06/10 19:03:13 nichols Exp $ Thomas Nichols, Marko Wilke
    if nargin < 4, u = .001; end
    if nargin < 5, alpha = .05; end
    if iscell(im), im = char(im); end

    assert(ismember(direct, {'+' '-' '+/-'}));

    V = spm_vol(im);
    Ymask = spm_read_vols(V);
    DIM = V.dim';
    [X,Y,Z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ0 = [X(:)';Y(:)';Z(:)'];

    %% Get Design Variable %%
    if exist([fileparts(im) filesep 'I.mat'],'file')
        matfile = [fileparts(im) filesep 'I.mat'];
        flexflag = 1;
    elseif exist([fileparts(im) filesep 'SPM.mat'],'file')
        matfile = [fileparts(im) filesep 'SPM.mat'];
        flexflag = 0;
    else
        disp('Could not find an SPM.mat or I.mat variable, exiting.'); extent = []; info = []; return
    end

    %% Load and Compute Params %%
    if flexflag % GLMFLEX

        II = load(matfile);
        try
            maskhdr = spm_vol(fullfile(fileparts(im), 'mask.nii'));
        catch
            maskhdr = spm_vol(fullfile(fileparts(im), 'mask.nii.gz'));
        end

        FWHM    = II.I.FWHM{1};
        R       = spm_resels_vol(maskhdr,FWHM)';
        VOX2RES = 1/prod(FWHM(~isinf(FWHM)));% voxels to resels

    else % SPM

        load(matfile);
        R    = SPM.xVol.R;
        VOX2RES  = 1/prod(SPM.xVol.FWHM(~isinf(SPM.xVol.FWHM))); %-voxels to resels

    end

    %% SPM METHOD (spm_uc_clusterFDR)
    thresh         = spm_invTcdf(1-u, df);
    di             = find(strcmpi({'+' '-' '+/-'}, direct));
    idx            = Ymask(:)~=0;
    Z              = Ymask(idx)';
    if di==2, Z = Z*-1; end
    XYZ            = XYZ0(:,idx);
    [fdrk,Pc,fwek] = spm_uc_clusterFDR(alpha, [1 df], 'T', R, 1, Z, XYZ, VOX2RES, thresh);
end



function [clustsize, clustidx]   = bspm_getclustidx(rawol, u, k)
    %
    % Get cluster indices after running bspm_cluster_correct(). Shamelessly stolen from bspmview
    %
    % INPUT:
    % im = raw volumes from the tmap .nii file, thresholded appropriately with threshold u
    % u = p-value threshold that was passed to bspm_cluster_correct(), e.g. 0.001
    % k = cluster extent threshold (FWE or FDR), result from bspm_cluster_correct()
    %
    % OUTPUT:
    % 

    % raw data to XYZ
    DIM         = size(rawol);
    [X,Y,Z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ         = [X(:)';Y(:)';Z(:)'];
    pos         = zeros(1, size(XYZ, 2));
    neg         = pos;
    clustidx    = zeros(3, size(XYZ, 2));

    % positive
    supra = (rawol(:)>u)';
    if sum(supra)
        tmp         = spm_clusters(XYZ(:, supra));
        clbin       = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
        pos(supra)  = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
        clustidx(1,supra) = tmp;
    end
    pos(pos < k)    = 0;

    % negative
    rawol = rawol*-1;
    supra = (rawol(:)>u)';
    if sum(supra)
        tmp      = spm_clusters(XYZ(:, supra));
        clbin      = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
        neg(supra) = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
        clustidx(2,supra) = tmp;
    end
    neg(neg < k) = 0;

    % both
    clustsize       = [pos; neg];
    clustsize(3,:)  = sum(clustsize);
    clustidx(3,:)   = sum(clustidx);
end




function PEAK = bspm_getmaxima(Z, XYZ, M, Dis, Num)
    [N,Z,XYZ,A,L]       = spm_max(Z,XYZ);
    XYZmm               = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
    npeak   = 0;
    PEAK    = [];
    while numel(find(isfinite(Z)))
        %-Find largest remaining local maximum
        %------------------------------------------------------------------
        [U,i]   = max(Z);            %-largest maxima
        j       = find(A == A(i));   %-maxima in cluster
        npeak   = npeak + 1;         %-number of peaks
        extent  = N(i);
        PEAK(npeak,:) = [extent U XYZmm(:,i)']; %-assign peak
        %-Print Num secondary maxima (> Dis mm apart)
        %------------------------------------------------------------------
        [l,q] = sort(-Z(j));                              % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = j(q(i));
            if min(sqrt(sum((XYZmm(:,D)-repmat(XYZmm(:,d),1,size(D,2))).^2)))>Dis
                if length(D) < Num
                    D          = [D d];
                    npeak   = npeak + 1;         %-number of peaks
                    PEAK(npeak,:) = [extent Z(d) XYZmm(:,d)']; %-assign peak
                end
            end
        end
        Z(j) = NaN;     % Set local maxima to NaN
    end
end



function [atlaslabels, atlas0] = bspm_setatlas(tmap_filename, atlas_dirpath, atlas_name)
    % setatlas from bspmview
    %
    atlas_vol = fullfile(atlas_dirpath, sprintf('%s_Atlas_Map.nii', atlas_name));
    if ~exist(atlas_vol, 'file')
        atlas_vol = fullfile(atlas_dirpath, sprintf('%s_Atlas_Map.nii.gz', atlas_name));
    end
    atlas_labels    = fullfile(atlas_dirpath, sprintf('%s_Atlas_Labels.mat', atlas_name));
    int             = 0; % 0=Nearest Neighbor, 1=Trilinear(default)
    atlasvol        = bspm_reslice_image(atlas_vol, tmap_filename, int);
    atlasvol        = single(round(atlasvol(:)))';
    load(atlas_labels);
    atlaslabels   = atlas;
    atlas0        = atlasvol;
end


function [out, outmat] = bspm_reslice_image(in, ref, int)
    % Most of the code is adapted from rest_Reslice in REST toolbox:
    % Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.
    % State Key Laboratory of Cognitive Neuroscience and Learning
    % Beijing Normal University, China, 100875
    % int:        interpolation, 0=Nearest Neighbor, 1=Trilinear(default)
    if nargin<3, int = 1; end
    if nargin<2, display('USAGE: [out, outmat] = reslice_image(infile, ref, SourceHead, int)'); return; end
    if iscell(ref), ref = char(ref); end
    if iscell(in), in = char(in); end
    % read in reference image
    RefHead = spm_vol(ref);
    mat=RefHead.mat;
    dim=RefHead.dim;
    SourceHead = spm_vol(in);
    [x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
    d       = [int*[1 1 1]' [1 1 0]'];
    C       = spm_bsplinc(SourceHead, d);
    v       = zeros(dim);
    M       = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
    y1      = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
    y2      = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
    y3      = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
    out     = spm_bsplins(C, y1,y2,y3, d);
    tiny = 5e-2; % From spm_vol_utils.c
    Mask = true(size(y1));
    Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
    Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
    Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
    out(~Mask) = 0;
    outmat = mat;
end



function [regionname, regionidx] = bspm_getregionnames(xyz, atlaslabels, atlas0, XYZmm0)
    if size(xyz,1)~=3, xyz = xyz'; end
    regionname  = repmat({'Location not in atlas'}, size(xyz, 2), 1);
    regionidx   = zeros(size(xyz, 2), 1);
    for i = 1:size(xyz,2)
        regionidx(i) = atlas0(spm_XYZreg('FindXYZ', xyz(:,i), XYZmm0));
        if regionidx(i)
            regionname{i} = atlaslabels.label{atlaslabels.id==regionidx(i)};
        end
    end
end


