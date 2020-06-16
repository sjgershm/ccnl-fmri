function [C, H, T, P, all_subject_coef, all_subject_z] = ccnl_match_rdms(Neural, Behavioral, control, corr_type)

    % Second-order RSA. Compute rank correlation between behavioral and neural RDMs, as well
    % as a group-level t-statistic showing whether that correlation is significant across subjects.
    %
    % USAGE:
    %   [C, H, T, P, all_subject_coef, all_subject_z] = ccnl_match_rdms(Neural, Behavioral, control, [corr_type])
    %
    % INPUT:
    %   Neural - [nROIs] struct array with neural RDMs, as output by ccnl_searchlight_rdms
    %   Behavioral - [nModels] struct array with neural RDMs, as output by ccnl_behavioral_rdms
    %   control - indices of models in Behavioral which to use as controls
    %   corr_type - (optional) type of (rank) correlation between RDMs (default: Spearman)
    %
    % OUTPUT:
    %   C - [nROIs x nModels] matrix of correlation coefficients (averaged across subjects)
    %   H - [nROIs x nModels] matrix of hypothesis outcomes
    %   T - [nROIs x nModels] matrix of t-statistics
    %   P - [nROIs x nModels] matrix of p-values
    %   all_subject_coef - [nROIs x nModels x nSubjects] matrix of (rank) correlations
    %   all_subject_z - [nROIs x nModels x nSubjects] matrix of (Fisher) z-transformed (rank) correlations (suitable for t-tests)
    %
    % Momchil Tomov, Oct 2018

    if ~exist('corr_type', 'var')
        corr_type = 'Spearman';
    end

    C = nan(length(Neural), length(Behavioral)); % average coef for each ROI for each model
    H = nan(length(Neural), length(Behavioral)); % result of hypothesis test for each ROI for each model -- is the correlation significant?
    T = nan(length(Neural), length(Behavioral)); % t-value of hypothesis test for each ROI for each model
    P = nan(length(Neural), length(Behavioral)); % p-value of hypothesis test for each ROI for each model
    all_subject_coef = nan(length(Neural), length(Behavioral), length(Behavioral(1).subj));
    all_subject_z = nan(length(Neural), length(Behavioral), length(Behavioral(1).subj));

    disp('matching RDMs...');
    tic

    for n = 1:length(Neural)
        fprintf('ROI #%d\n', n);

        models_subjs_coefs = []; % Spearman coefs: row = model, col = subject
        for m = 1:length(Behavioral)
            
            % Compute a Spearman's rank correlation for each subject separately
            %
            subjs_coefs = []; % Spearman coefs: col = subject, one coef per
            for s = 1:length(Behavioral(1).subj)
                assert(isequal(size(Neural(n).subj(s).RDM), size(Behavioral(m).subj(s).RDM)), 'Neural and behavioral RDMs should be equal -- check rsa.which_trials');

                upper = logical(triu(ones(size(Neural(n).subj(s).RDM)), 1)); % only take values above main diagonal
                upper = upper & Behavioral(m).subj(s).partition_RDM; % don't compare within the same partition (run)! BOLD is VERY autocorrelated

                neural_RDM = Neural(n).subj(s).RDM(upper); 
                behavioral_RDM = Behavioral(m).subj(s).RDM(upper);

                % control RDMs
                control_RDMs = [];
                if ~ismember(m, control) % do not correct control RDMs for themselves (you get nothing)
                    for i = 1:length(control)
                        control_RDM = Behavioral(control(i)).subj(s).RDM(upper);
                        control_RDMs(:, i) = control_RDM;
                    end
                end

                % rank correlate RDMs
                switch lower(corr_type)
                    case {'spearman', 'pearson'}
                        % one of built-ins

                        if ~isempty(control_RDMs)
                            % partial correlation
                            subj_coef = partialcorr(neural_RDM, behavioral_RDM, control_RDMs, 'type', corr_type);
                        else
                            subj_coef = corr(neural_RDM, behavioral_RDM, 'type', corr_type);
                        end

                        assert(~isequal(lower(corr_type), 'kendall'), 'Kendall standardization is not atanh');
                        subj_z = atanh(subj_coef); % fisher z-transform

                    case {'ktaub'}

                        % Kendall's tau b, adjusting for ties, unlike Spearman or MATLAB's Kendall

                        assert(isempty(control_RDMs), 'partial correlation not supported for ktaub yet');
                        [subj_coef, ~, ~, ~, subj_z] = ktaub([neural_RDM behavioral_RDM], 0.05);

                    otherwise
                        assert(false, 'invalid corr_type');
                end
                all_subject_coef(n, m, s) = subj_coef;
                all_subject_z(n, m, s) = subj_z;

            end

        end

        % Group-level analysis -- collapse across subjects
        % for each model, do a t-test of the correlation coefficients between
        % each subject's neural RDM and the corresponding model RDM
        %
        [h, ps, ci, stats] = ttest(permute(all_subject_z(n, :, :), [3 2 1])); % t-test across subjects

        C(n,:) = mean(permute(all_subject_coef(n, :, :), [3 2 1]), 1); % average across subjects
        H(n,:) = h;
        P(n,:) = ps;
        T(n,:) = stats.tstat;
    end

    toc

