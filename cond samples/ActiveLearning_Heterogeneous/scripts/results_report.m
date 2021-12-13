%% batch script: data prep and results reporting for Active Learning Heterogeneous Sample
%a_set_up_paths_always_run_first
b_reload_saved_matfiles

%% data prep (uncomment if needs to be redone)
%z_batch_publish_image_prep_and_qc

%results reporting
%% SVM results
% had to go into plugin_svm_contrasts_get_results_per_subject and 
% uncomment line216/ comment line217 to make this work. don't really know if result is valid
% did not have to do this for 0188 sample, why? different N but design is the same??
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