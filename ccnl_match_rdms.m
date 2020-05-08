function [Rho, H, T, P, all_subject_rhos] = ccnl_match_rdms(Neural, Behavioral, control)

    % Second-order RSA. Compute rank correlation between behavioral and neural RDMs, as well
    % as a group-level t-statistic showing whether that correlation is significant across subjects.
    %
    % USAGE:
    %   [Rho, H, T, P, all_subject_rhos] = ccnl_match_rdms(Neural, Behavioral, control)
    %
    % INPUT:
    %   Neural - [nROIs] struct array with neural RDMs, as output by ccnl_searchlight_rdms
    %   Behavioral - [nModels] struct array with neural RDMs, as output by ccnl_behavioral_rdms
    %   control - indices of models in Behavioral which to use as controls
    %
    % OUTPUT:
    %   Rho - [nROIs x nModels] matrix of Spearman rank correlations (averaged across subjects)
    %   H - [nROIs x nModels] matrix of hypothesis outcomes
    %   T - [nROIs x nModels] matrix of t-statistics
    %   P - [nROIs x nModels] matrix of p-values
    %   all_subject_rhos - [nROIs x nModels x nSubjects] matrix of Spearman rank correlations
    %
    % Momchil Tomov, Oct 2018


    Rho = []; % average Spearman's rho for each ROI for each model
    H = []; % result of hypothesis test for each ROI for each model -- is the correlation significant?
    T = []; % t-value of hypothesis test for each ROI for each model
    P = []; % p-value of hypothesis test for each ROI for each model

    for n = 1:length(Neural)
        fprintf('ROI #%d\n', n);

        models_subjs_rhos = []; % Spearman rhos: row = model, col = subject
        for m = 1:length(Behavioral)
            
            % Compute a Spearman's rank correlation for each subject separately
            %
            subjs_rhos = []; % Spearman rhos: col = subject, one rho per
            for s = 1:length(Behavioral(1).subj)
                assert(isequal(size(Neural(n).subj(s).RDM), size(Behavioral(m).subj(s).RDM)), 'Neural and behavioral RDMs should be equal -- check rsa.which_trials');

                upper = logical(triu(ones(size(Neural(n).subj(s).RDM)), 1)); % only take values above main diagonal
                upper = upper & Behavioral(m).subj(s).run_RDM; % don't compare within the same run! BOLD is VERY autocorrelated

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
                if ~isempty(control_RDMs)
                    subj_rho = partialcorr(neural_RDM, behavioral_RDM, control_RDMs, 'type', 'Spearman');
                else
                    subj_rho = corr(neural_RDM, behavioral_RDM, 'type', 'Spearman');
                end
                subjs_rhos = [subjs_rhos, subj_rho];
                all_subject_rhos(n, m, s) = subj_rho;            

            end

            models_subjs_rhos = [models_subjs_rhos; subjs_rhos]; % for group-level analysis
        end

        % Group-level analysis -- collapse across subjects
        % for each model, do a t-test of the correlation coefficients between
        % each subject's neural RDM and the corresponding model RDM
        %
        fisher_models_rhos = atanh(models_subjs_rhos);
        [h, ps, ci, stats] = ttest(fisher_models_rhos');

        Rho(n,:) = mean(fisher_models_rhos');
        H(n,:) = h;
        P(n,:) = ps;
        T(n,:) = stats.tstat;
    end

