%% batch script: data prep and results reporting for Active Learning Homogeneous Sample

%% data prep (uncomment if needs to be redone)

% a_set_up_paths_always_run_first
% z_batch_publish_image_prep_and_qc

%results reporting
%% SVM results
c2_SVM_contrasts_masked                   
% Show SVM results (from prep_3b_run_SVMs_on_contrasts_and_save)

%% Univariate results
c_univariate_contrast_maps                        % Show voxel-wise maps for all contrasts, FDR-corrected and uncorrected; tables
c3_univariate_contrast_maps_scaleddata            % Show voxel-wise contrast maps for scaled/cleaned data
c4_univ_contrast_wedge_plot_and_lateralization_test   % Test contrasts and lateralization across 16 large-scale rsFMRI networks with L/R divisions


%% Compare to emotion signatures
d5_emotion_signature_responses
d12_kragel_emotion_signature_responses
% d13_kragel_emotion_riverplots
d14_kragel_emotion_signature_similarity_barplots

%% Compare to Buckner networks
f2_bucknerlab_network_barplots                    % Contrasts analyzed according to 7 Buckner Lab cortical networks from Yeo et al. 2011
f2_bucknerlab_network_wedgeplots
% g2_bucknerlab_network_riverplots

%% Compare to empathy
d6_empathy_signature_responses

%%