function table = ccnl_results_table(varargin)
    %
    % Same as bspmview -> Show Results Table (after cluster FWE correction). Except useful.
    % Given a contrast, extract all the activation clusters from the t-map after cluster FWE correction.
    % Returns result as a table with the cluster names, extents, t-statistics, etc.
    % As a bonus, prints out the table in LaTeX format.
    %
    % USAGE:
    %   table = ccnl_results_table(atlas_name, method, EXPT, model, contrast, p, direct, alpha, Dis, Num)
    %
    % EXAMPLES:
    %   table = ccnl_results_table('AAL2', 'peak', exploration_expt(), 1, 'RR - SS')
    %   table = ccnl_results_table('Brodmann', 'vote', exploration_expt(), 1, 'RR - SS', 0.001, '+', 0.05)
    %
    % INPUT:
    %   atlas_name = name of atlas to use; options are 'AAL2', 'AnatomyToolbox', 'HarvardOxford-maxprob-thr0', 'Brodmann', 'Talairach' (defaults to 'AAL2')
    %   method = how to pick the cluster names, should be one of 'peak', 'vote', or 'all' (bspmview default is 'peak')
    %   ...the rest of the parameters are the same as in ccnl_extract_clusters()
    %
    % OUTPUT:
    %   table = the bspmview results table
    % 

    [curdir, ~, ~] = fileparts(mfilename('fullpath')); 
    atlas_dirpath = fullfile(curdir, 'atlases');

    atlas_name = varargin{1};
    method = varargin{2};
    assert(ismember(atlas_name, {'AnatomyToolbox', 'AAL2', 'HarvardOxford-maxprob-thr0', 'Talairach', 'Brodmann'}));
    assert(ismember(method, {'peak', 'vote', 'all'}));

    % extract the clusters
    %
    [V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(varargin{3:end});

    % load atlas
    %
    [atlaslabels, atlas] = bspm_setatlas(spmT, atlas_dirpath, atlas_name);
    atlas = reshape(atlas, size(Y)); % voxel --> num
    map = containers.Map(atlaslabels.id, atlaslabels.label); % num --> label

    % load BA atlas separately
    %
    [~, BAatlas] = bspm_setatlas(spmT, atlas_dirpath, 'Brodmann');
    BAatlas = reshape(BAatlas, size(Y)); % voxel --> BA

    % label the clusters
    %
    new_region = {};

    for i = 1:size(region, 1) 
        % get peak voxel
        x = cor(i,1);
        y = cor(i,2);
        z = cor(i,3);
        assert(immse(stat(i), Y(x,y,z)) < 1e-6);
        
        % set brodmann area
        if BAatlas(x, y, z)
            BA = num2str(BAatlas(x, y, z));
        else
            BA = '';
        end

        % get cluser mask and voxels
        clust_idx = CI(x,y,z);
        mask = CI == clust_idx;
        voxels = find(mask);
        if ~isequal(method, 'peak')
            assert(numel(voxels) == extent(i), 'Please set direct to ''+'' or ''-''; sometimes it does not work with ''+/-'''); % <-- doesn't work with +/- ; do one at a time
        end

        % see which region each voxel "votes" for
        %
        [x, y, z] = ind2sub(size(mask), voxels);
        votes = zeros(1, max(atlas(:)));
        sign_votes = zeros(1, max(atlas(:))); % whether it's L or R on average
        for j = 1:numel(x)
            if atlas(x(j),y(j),z(j)) > 0
                idx = atlas(x(j),y(j),z(j));
                votes(idx) = votes(idx) + 1;
                
                m = cor2mni([x(j) y(j) z(j)], V.mat);
                if m(1) < 0
                    sign_votes(idx) = sign_votes(idx) - 1; % L voxel
                else
                    sign_votes(idx) = sign_votes(idx) + 1; % R voxel
                end
            end
        end

        % get the cluster label(s)
        %
        switch method
            case 'peak'
                % just use the peak voxel
                %
                if isKey(map, atlas(cor(i,1),cor(i,2),cor(i,3)))
                    label = map(atlas(cor(i,1),cor(i,2),cor(i,3)));
                    new_region{i} = atlas_label_to_roi_name(atlas_name, label, mni(i,:)); % convert to nice anatomical name
                else
                    new_region{i} = '';
                end
                
            case 'vote'
                % get the most popular anatomical region in the cluster
                % each voxel gets one vote
                %
                if sum(votes) == 0
                    new_region{i} = ''; % no region :(
                else
                    [~, idx] = max(votes);
                    label = map(idx);
                    new_region{i} = atlas_label_to_roi_name(atlas_name, label, [sign_votes(idx) 0 0]); % convert to nice anatomical name
                end
                
            case 'all'
                % List the regions of all voxels
                %
                if sum(votes) == 0
                    new_region{i} = ''; % no region :(
                else
                    [~, idxs] = find(votes);
                    [~, sorted_idxs] = sort(votes(idxs), 'descend');

                    rois = {};
                    for j = 1:numel(sorted_idxs)
                        idx = idxs(sorted_idxs(j));
                        assert(votes(idx) > 0);

                        label = map(idx);
                        roi = atlas_label_to_roi_name(atlas_name, label, [sign_votes(idx) 0 0]);
                        roi = sprintf('%s (%.2f\\%%)', roi, 100 * votes(idx) / sum(votes));
                        rois = [rois; {roi}];
                    end
                    new_region{i} = rois;
                end
                
            otherwise
                assert(false, 'method should be one of these');
        end
        
        % print for LaTeX
        %
        sign = '';
        if i == 1 || (stat(i) < 0) ~= (stat(i - 1) < 0) % sign changed
            if stat(i) < 0
                sign = 'Negative';
            else
                sign = 'Positive';
            end
        end
        
        region_latex = new_region{i};
        if iscell(new_region{i}) % multiple regions in this cluster
            region_latex = strjoin(new_region{i}, ' \\\\\n & ');
        else
            region_latex = new_region{i};
        end
        fprintf('%s & %s & %s & %d & %.3f & %d %d %d \\\\\n', sign, region_latex, BA, extent(i), stat(i), mni(i,1), mni(i,2), mni(i,3));
        if isequal(method, 'all')
            fprintf('\\hline\n');
        end
        
        % output table
        %
        table{i, 1} = sign;
        table{i, 2} = new_region{i};
        table{i, 3} = BA;
        table{i, 4} = extent(i);
        table{i, 5} = stat(i);
        table{i, 6} = mni(i,1);
        table{i, 7} = mni(i,2);
        table{i, 8} = mni(i,3);

        
    end


    %new_region'
end


function roi = atlas_label_to_roi_name(atlas_name, roi_label, mni)

    % Function that converts AAL2 labels from e.g. bspmview tables into real anatomical ROI names.
    % 
    % INPUT:
    % atlas_name = name of the atlas to use
    % roi_label = ROI label from that atlas, e.g. 'Angular_R' for AAL2
    % mni = MNI coordinates to infer laterality
    %
    % OUTPUT:
    % roi = the actual anatomical name, e.g. 'Angular gyrus (R)'. If not round,
    %       returns the passed roi_label
    %

    roi = roi_label;

    switch atlas_name
        case 'AAL2'
            roi = aal2_label_to_roi_name(roi, mni);

        case 'AnatomyToolbox'
            hemi = '';
            if startsWith(roi, 'L ')
                roi = roi(3:end);
                hemi = 'L';
            elseif startsWith(roi, 'R ')
                roi = roi(3:end);
                hemi = 'R';
            end

            space = find(roi == ' ' | roi == '-');
            if ~isempty(space)
                roi = [roi(1:space), lower(roi(space+1:end))];
            end

            if ~isempty(hemi)
                roi = [roi, ' (', hemi, ')'];
            end

        case 'HarvardOxford-maxprob-thr0' 
            hemi = '';
            if  mni(1) < 0
                roi = [roi, ' (L)'];
            else
                roi = [roi, ' (R)'];
            end
            
        case 'Talairach'
            roi = roi;

        case 'Brodmann' 
            hemi = '';
            if  mni(1) < 0
                roi = [roi, ' (L)'];
            else
                roi = [roi, ' (R)'];
            end
            
        otherwise
            assert(false, 'atlas_name be one of the above');
    end
end


function roi = aal2_label_to_roi_name(roi_label, mni)

    % Function that converts AAL2 labels from e.g. bspmview tables into real anatomical ROI names.
    % 
    % INPUT:
    % roi_label = ROI label from AAL2 atlas as output by bspmview, e.g. 'Angular_R'.
    % mni = optional MNI coordinates to include laterality
    %
    % OUTPUT:
    % roi = the actual anatomical name, e.g. 'Angular gyrus' or 'Angular gyrus
    %       (R)'. If not round, returns the passed roi_label
    %

    % Anatomical region name, AAL2 label
    % Table 2 from Rolls et al., Implementation of a new parcellation of the orbitofrontal cortex in the
    % automated anatomical labeling atlas (NeuroImage, 2015)
    %
    aal2_labels = {...
    'Precentral gyrus',  'Precentral' ;
    'Postcentral gyrus',  'Postcentral' ;
    'Rolandic operculum',  'Rolandic_Oper' ;
    'Superior frontal gyrus, dorsolateral',  'Frontal_Sup' ;
    'Middle frontal gyrus',  'Frontal_Mid' ;
    'IFG pars opercularis',  'Frontal_Inf_Oper' ;
    'IFG pars triangularis',  'Frontal_Inf_Tri' ;
    'Superior frontal gyrus, medial',  'Frontal_Sup_Med' ;
    'Supplementary motor area',  'Supp_Motor_Area' ;
    'Paracentral lobule',  'Paracentral_Lobule' ;
    'Superior frontal gyrus, medial orbital',  'Frontal_Med_Orb' ;
    'IFG pars orbitalis',  'Frontal_Inf_Orb' ;
    'Gyrus rectus',  'Rectus' ;
    'Medial orbital gyrus',  'OFCmed' ;
    'Anterior orbital gyrus',  'OFCant' ;
    'Posterior orbital gyrus',  'OFCpost' ;
    'Lateral orbital gyrus',  'OFClat' ;
    'Olfactory cortex',  'Olfactory' ;
    'Superior temporal gyrus',  'Temporal_Sup' ;
    'Heschl''s gyrus',  'Heschl' ;
    'Middle temporal gyrus',  'Temporal_Mid' ;
    'Inferior temporal gyrus',  'Temporal_Inf' ;
    'Superior parietal gyrus',  'Parietal_Sup' ;
    'Inferior parietal gyrus, excluding supramarginal and angular gyri',  'Parietal_Inf' ;
    'Angular gyrus',  'Angular' ;
    'Supramarginal gyrus',  'SupraMarginal',;
    'Precuneus',  'Precuneus' ;
    'Superior occipital gyrus',  'Occipital_Sup' ;
    'Middle occipital gyrus',  'Occipital_Mid' ;
    'Inferior occipital gyrus',  'Occipital_Inf' ;
    'Cuneus',  'Cuneus' ;
    'Calcarine fissure and surrounding cortex',  'Calcarine' ;
    'Lingual gyrus',  'Lingual' ;
    'Fusiform gyrus',  'Fusiform' ;
    'Temporal pole: superior temporal gyrus',  'Temporal_Pole_Sup' ;
    'Temporal pole: middle temporal gyrus',  'Temporal_Pole_Mid' ;
    'Anterior cingulate \& paracingulate gyri',  'Cingulate_Ant' ;
    'Middle cingulate \& paracingulate gyri',  'Cingulate_Mid' ;
    'Posterior cingulate gyrus',  'Cingulate_Post' ;
    'Hippocampus',  'Hippocampus';
    'Parahippocampal gyrus',  'ParaHippocampal' ;
    'Insula',  'Insula' ;
    'Amygdala',  'Amygdala' ;
    'Caudate nucleus',  'Caudate' ;
    'Lenticular nucleus, Putamen',  'Putamen' ;
    'Lenticular nucleus, Pallidum',  'Pallidum' ;
    'Thalamus',  'Thalamus';
    'Cerebellum', 'Cerebelum';
    'Vermis', 'Vermis'};

    roi = roi_label;

    for j=1:size(aal2_labels, 1)
        if startsWith(roi_label, aal2_labels{j, 2})
            roi = aal2_labels{j, 1};
            if exist('mni', 'var') % optionally include laterality
                if mni(1) < 0
                    roi = [roi, ' (L)'];
                else
                    roi = [roi, ' (R)'];
                end
            end
            break;
        end
    end
end



%
% TODO dedupe with ccnl_extract_clusters
%

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
