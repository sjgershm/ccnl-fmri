function [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx, subjects)

    % Computes behavioral RDMs.
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx [, subjects])
    %
    % EXAMPLE:
    %   [Behavioral, control] = ccnl_behavioral_rdms(exploration_expt(), 1);
    %   showRDMs(Behavioral(1).subj);
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   subjects - (optional) list of subjects
    %
    % OUTPUTS:
    %   Behavioral - struct array, one element per model, with the fields:
    %      .name - model name (e.g. 'condition')
    %      .D - dimension of features
    %      .is_control - boolean whether this is a control model
    %      .distance_measure - name (e.g. 'cosine') or function handler to be used as a distance measure for the RDMs (passed to pdist, see MATLAB documentation)
    %      .subj - struct array with subject RDMs for given model; has following fields:
    %         .features - [nTrials x D] feature vector
    %         .runs - [nTrials x 1] run id vector
    %         .RDM - [nTrials x nTrials] RDM based on features
    %         .run_RDM - [nTrials x nTrials] run logical RDM (to compare across runs only)
    %   control - indices of models in Behavioral that are controls
    %
    % Momchil Tomov, Sep 2018

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end
    
    % for each subject
    for s = 1:length(subjects)
        subj = subjects(s);

        rsa = EXPT.create_rsa(rsa_idx, subj);
        n = length(rsa.model);

        % initialize model metadata for all models
        for i = 1:n % for each RSA model
            if s == 1
                Behavioral(i).D = size(rsa.model(i).features, 2);
                Behavioral(i).name = rsa.model(i).name;
                Behavioral(i).is_control = rsa.model(i).is_control;
                Behavioral(i).distance_measure = rsa.model(i).distance_measure;
            else
                assert(Behavioral(i).D == size(rsa.model(i).features, 2));
                assert(isequal(Behavioral(i).name, rsa.model(i).name));
                assert(Behavioral(i).is_control == rsa.model(i).is_control);
            end
            Behavioral(i).subj(s).name = Behavioral(i).name;
            Behavioral(i).subj(s).features = rsa.model(i).features;
            Behavioral(i).subj(s).runs = rsa.model(i).runs;
            Behavioral(i).subj(s).id = subj;
        end
    end

    % compute RDMs
    for i = 1:length(Behavioral)
        for s = 1:length(Behavioral(i).subj)
            Behavioral(i).subj(s).RDM = squareRDMs(pdist(Behavioral(i).subj(s).features, Behavioral(i).distance_measure));
            assert(sum(any(isnan(Behavioral(i).subj(s).RDM))) == 0, 'Found NaNs in RDM -- should never happen');

            Behavioral(i).subj(s).run_RDM = logical(squareRDMs(pdist(Behavioral(i).subj(s).runs, @(r1, r2) r1 ~= r2)));
            assert(sum(any(isnan(Behavioral(i).subj(s).run_RDM))) == 0, 'Found NaNs in run RDM -- should never happen');
        end
    end

    control = find([Behavioral.is_control]);
end

