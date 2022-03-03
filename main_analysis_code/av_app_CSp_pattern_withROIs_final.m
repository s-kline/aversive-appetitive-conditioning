%% author: s-kline
% analysis for paper comparing appetitive and aversive conditioning
clear all, close all
%% -- prep -- 
%% set up directories
base_dir = 'D:\Dropbox\aversive appetitive conditioning\github_repo';
% multivariate pattern dirs 
meta_dir = fullfile(base_dir, '\n_signature_fear_conditioning\');
% appetitive cond sample dirs
MID_0188_dir = fullfile(base_dir,'cond samples\ActiveLearning_Homogeneous\'); % active hom sample
MID_0205_dir = fullfile(base_dir,'cond samples\ActiveLearning_Heterogeneous\'); % active het sample
Acq_0187_dir = fullfile(base_dir,'cond samples\PassiveLearning_Heterogeneous\'); % passive het sample

%% load data
% -- neural fear conditiong signature map 
% meta analysis, Fullana et al 2016 
meta_fearcond_obj = statistic_image(fullfile(meta_dir, 'MyMean_z.nii.gz')); % map of activations and deactivations related to aversive CS+
color_meta_fearcond = [224  0  97]/255; % ruby (E00061)

% -------------------------------------------------------------------------
% -- appetitive conditioning data
% Active Learning/Homogeneous Sample (0188ok study, monetary incentive delay)
% Kruse et al 2018
color_CSpGTCSm_0188 = [76 201 240]/255; % vivid sky blue (4CC9F0)
MID_0188_struct = importdata(fullfile(MID_0188_dir, '\results\image_names_and_setup.mat'));
MID_0188_DATA_OBJ = importdata(fullfile(MID_0188_dir, '\results\data_objects.mat'));
MID_0188_CONTRASTS_OBJ = importdata(fullfile(MID_0188_dir, '\results\contrast_data_objects.mat'));
MID_0188_svm_stats_results = importdata(fullfile(MID_0188_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
MID_0188_DAT = MID_0188_struct.DAT;
CSplus_0188 = fmri_data(MID_0188_DATA_OBJ{1, 3});
CSminus_0188 = fmri_data(MID_0188_DATA_OBJ{1, 6});
CSpGTCSm_0188_con_obj = MID_0188_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1};
CSpGTCSm_0188_svm_obj = MID_0188_svm_stats_results{1,1}.weight_obj;
% -------------------------------------------------------------------------
% Active Learning/Heterogeneous Sample (0205ok study, monetary incentive delay)
% Kruse et al 2020
color_CSpGTCSm_0205 = [58 12 163]/255; % blue (trypan blue 3A0CA3)
MID_0205_struct = importdata(fullfile(MID_0205_dir, '\results\image_names_and_setup.mat'));
MID_0205_DATA_OBJ = importdata(fullfile(MID_0205_dir, '\results\data_objects.mat'));
MID_0205_CONTRASTS_OBJ = importdata(fullfile(MID_0205_dir, '\results\contrast_data_objects.mat'));
MID_0205_svm_stats_results = importdata(fullfile(MID_0205_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
MID_0205_DAT = MID_0205_struct.DAT;
CSplus_0205 = fmri_data(MID_0205_DATA_OBJ{1, 3});
CSminus_0205 = fmri_data(MID_0205_DATA_OBJ{1, 6});
CSpGTCSm_0205_con_obj = MID_0205_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1};
CSpGTCSm_0205_svm_obj = MID_0205_svm_stats_results{1,1}.weight_obj;
% -------------------------------------------------------------------------
% Passive Learning/Heterogeneous Sample (0187tk study, passive task)
% Tapia León et al 2019
color_CSpGTCSm_0187 = [247 37 133]/255; % flickr pink (F72585)
Acq_0187_struct = importdata(fullfile(Acq_0187_dir, '\results\image_names_and_setup.mat'));
Acq_0187_DATA_OBJ = importdata(fullfile(Acq_0187_dir, '\results\data_objects.mat'));
Acq_0187_CONTRASTS_OBJ = importdata(fullfile(Acq_0187_dir, '\results\contrast_data_objects.mat'));
Acq_0187_svm_stats_results = importdata(fullfile(Acq_0187_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
Acq_0187_DAT = Acq_0187_struct.DAT;
CSplus_0187 = fmri_data(Acq_0187_DATA_OBJ{1, 1});
CSminus_0187 = fmri_data(Acq_0187_DATA_OBJ{1, 2});
CSpGTCSm_0187_con_obj = Acq_0187_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1};
CSpGTCSm_0187_svm_obj = Acq_0187_svm_stats_results{1,1}.weight_obj;
% -------------------------------------------------------------------------

%% prepare masks for ROI analysis
% load atlas & regions of interest
atlas_obj = load_atlas('canlab2018_2mm');
basal_ganglia = select_atlas_subset(atlas_obj, 'labels_2', {'Basal_ganglia'});
% Nucleus Accumbens
nacc = select_atlas_subset(basal_ganglia, 'labels', {'NAC_L', 'NAC_R'}, 'flatten');
nacc = resample_space(nacc, CSpGTCSm_0188_con_obj); % mask space needs to be resmpled to data space (the samples all have the same space)
nacc.fullpath='nacc.nii';nacc.write('overwrite');
% Caudate Nucleus 
caudate = select_atlas_subset(basal_ganglia, 'labels', {'Caudate_Cp_L', 'Caudate_Cp_R', 'Caudate_Ca_L', 'Caudate_Ca_R'}, 'flatten');
caudate = resample_space(caudate, CSpGTCSm_0188_con_obj);
caudate.fullpath='caudate.nii';caudate.write('overwrite');
% Amygdala
amygdala = select_atlas_subset(atlas_obj, 'labels_2', {'Amygdala'}, 'flatten');
amygdala = resample_space(amygdala, CSpGTCSm_0188_con_obj);
amygdala.fullpath='amygdala.nii';amygdala.write('overwrite');
% Cerebellum
cerebellum = select_atlas_subset(atlas_obj, 'labels_2', {'Cerebellum'}, 'flatten');
cerebellum = resample_space(cerebellum, CSpGTCSm_0188_con_obj);
cerebellum.fullpath='cerebellum.nii';cerebellum.write('overwrite');
% -------------------------------------------------------------------------

%% mask pattern (only for figures, not analysis)
roi_meta_fearcond_obj_nacc = apply_mask(meta_fearcond_obj, nacc);
roi_meta_fearcond_obj_caudate = apply_mask(meta_fearcond_obj, caudate);
roi_meta_fearcond_obj_amygdala = apply_mask(meta_fearcond_obj, amygdala);
roi_meta_fearcond_obj_cerebellum = apply_mask(meta_fearcond_obj, cerebellum);

%% mask sample data
% Active Learning/Homogeneous Sample 
roi_CSpGTCSm_0188_con_obj_nacc = apply_mask(CSpGTCSm_0188_con_obj, nacc);
roi_CSpGTCSm_0188_con_obj_caudate = apply_mask(CSpGTCSm_0188_con_obj, caudate);
roi_CSpGTCSm_0188_con_obj_amygdala = apply_mask(CSpGTCSm_0188_con_obj, amygdala);
roi_CSpGTCSm_0188_con_obj_cerebellum = apply_mask(CSpGTCSm_0188_con_obj, cerebellum);
roi_CSplus_0188_nacc = apply_mask(CSplus_0188, nacc);
roi_CSminus_0188_nacc = apply_mask(CSminus_0188, nacc);
roi_CSplus_0188_caudate = apply_mask(CSplus_0188, caudate);
roi_CSminus_0188_caudate = apply_mask(CSminus_0188, caudate);
roi_CSplus_0188_amygdala = apply_mask(CSplus_0188, amygdala);
roi_CSminus_0188_amygdala = apply_mask(CSminus_0188, amygdala);
roi_CSplus_0188_cerebellum = apply_mask(CSplus_0188, cerebellum);
roi_CSminus_0188_cerebellum = apply_mask(CSminus_0188, cerebellum);

% -------------------------------------------------------------------------
% Active Learning/Heterogeneous Sample
roi_CSpGTCSm_0205_con_obj_nacc = apply_mask(CSpGTCSm_0205_con_obj, nacc);
roi_CSpGTCSm_0205_con_obj_caudate = apply_mask(CSpGTCSm_0205_con_obj, caudate);
roi_CSpGTCSm_0205_con_obj_amygdala = apply_mask(CSpGTCSm_0205_con_obj, amygdala);
roi_CSpGTCSm_0205_con_obj_cerebellum = apply_mask(CSpGTCSm_0205_con_obj, cerebellum);
roi_CSplus_0205_nacc = apply_mask(CSplus_0205, nacc);
roi_CSminus_0205_nacc = apply_mask(CSminus_0205, nacc);
roi_CSplus_0205_caudate = apply_mask(CSplus_0205, caudate);
roi_CSminus_0205_caudate = apply_mask(CSminus_0205, caudate);
roi_CSplus_0205_amygdala = apply_mask(CSplus_0205, amygdala);
roi_CSminus_0205_amygdala = apply_mask(CSminus_0205, amygdala);
roi_CSplus_0205_cerebellum = apply_mask(CSplus_0205, cerebellum);
roi_CSminus_0205_cerebellum = apply_mask(CSminus_0205, cerebellum);

% -------------------------------------------------------------------------
% Passive Learning/Heterogeneous Sample
roi_CSpGTCSm_0187_con_obj_nacc = apply_mask(CSpGTCSm_0187_con_obj, nacc);
roi_CSpGTCSm_0187_con_obj_caudate = apply_mask(CSpGTCSm_0187_con_obj, caudate);
roi_CSpGTCSm_0187_con_obj_amygdala = apply_mask(CSpGTCSm_0187_con_obj, amygdala);
roi_CSpGTCSm_0187_con_obj_cerebellum = apply_mask(CSpGTCSm_0187_con_obj, cerebellum);
roi_CSplus_0187_nacc = apply_mask(CSplus_0187, nacc);
roi_CSminus_0187_nacc = apply_mask(CSminus_0187, nacc);
roi_CSplus_0187_caudate = apply_mask(CSplus_0187, caudate);
roi_CSminus_0187_caudate = apply_mask(CSminus_0187, caudate);
roi_CSplus_0187_amygdala = apply_mask(CSplus_0187, amygdala);
roi_CSminus_0187_amygdala = apply_mask(CSminus_0187, amygdala);
roi_CSplus_0187_cerebellum = apply_mask(CSplus_0187, cerebellum);
roi_CSminus_0187_cerebellum = apply_mask(CSminus_0187, cerebellum);

% -------------------------------------------------------------------------
%% -- Main Analysis -- 
%% test meta fear conditioning signature on appetitive cond data -- 
% ROI analysis: masked data with ROIs and calculated similarity.
% the canlab_pattern_similarity doc implies, that the pattern should not be
% masked, because of how it deals with zeroes in statistic images vs data
% "Thus, this function treats values of 0 in DATA images, not pattern masks,
% as missing values, and excludes these voxels from both image and
% mask when calculating similarity."
metric = 'cosine_similarity';
%% Active Learning/Homogeneous Sample, 0188tk 
% similarity of CSpGTCSm contrast data and fear conditioning signature pattern
appCSpGTCSm0188_meta_fear_sim = apply_mask(CSpGTCSm_0188_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
CSplus_0188_meta_fear_expression = apply_mask(CSplus_0188,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
CSminus_0188_meta_fear_expression = apply_mask(CSminus_0188,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);

% similarity of data masked with ROIs with pattern
nacc_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_nacc,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
caudate_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_caudate,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
amygdala_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_amygdala,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
cerebellum_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_cerebellum,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
nacc_CSplus_0188_meta_fear_expression = apply_mask(roi_CSplus_0188_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
nacc_CSminus_0188_meta_fear_expression = apply_mask(roi_CSminus_0188_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSplus_0188_meta_fear_expression = apply_mask(roi_CSplus_0188_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSminus_0188_meta_fear_expression = apply_mask(roi_CSminus_0188_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSplus_0188_meta_fear_expression = apply_mask(roi_CSplus_0188_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSminus_0188_meta_fear_expression = apply_mask(roi_CSminus_0188_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSplus_0188_meta_fear_expression = apply_mask(roi_CSplus_0188_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSminus_0188_meta_fear_expression = apply_mask(roi_CSminus_0188_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);

% negative control patterns
appCSpGTCSm0188_cog_sim = apply_all_signatures(CSpGTCSm_0188_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pain_cog_emo');
                                           
appCSpGTCSm0188_emotion_sim = apply_all_signatures(CSpGTCSm_0188_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'kragelemotion');
                                           
appCSpGTCSm0188_pines_sim = apply_all_signatures(CSpGTCSm_0188_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pines');
                                           
appCSpGTCSm0188_stroop_sim = apply_all_signatures(CSpGTCSm_0188_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'stroop');
%% Active Learning/Heterogenoeus Sample, 0205ok 
% similarity of CSpGTCSm contrast data and fear conditioning signature pattern
appCSpGTCSm0205_meta_fear_sim = apply_mask(CSpGTCSm_0205_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
CSplus_0205_meta_fear_expression = apply_mask(CSplus_0205,replace_empty(meta_fearcond_obj),  'pattern_expression', metric);
CSminus_0205_meta_fear_expression = apply_mask(CSminus_0205,replace_empty(meta_fearcond_obj),  'pattern_expression', metric);

% similarity of data masked with ROIs with pattern
nacc_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_nacc,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
caudate_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_caudate,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
amygdala_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_amygdala,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
cerebellum_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_cerebellum,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
nacc_CSplus_0205_meta_fear_expression = apply_mask(roi_CSplus_0205_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
nacc_CSminus_0205_meta_fear_expression = apply_mask(roi_CSminus_0205_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSplus_0205_meta_fear_expression = apply_mask(roi_CSplus_0205_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSminus_0205_meta_fear_expression = apply_mask(roi_CSminus_0205_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSplus_0205_meta_fear_expression = apply_mask(roi_CSplus_0205_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSminus_0205_meta_fear_expression = apply_mask(roi_CSminus_0205_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSplus_0205_meta_fear_expression = apply_mask(roi_CSplus_0205_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSminus_0205_meta_fear_expression = apply_mask(roi_CSminus_0205_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);

% negative control patterns
appCSpGTCSm0205_cog_sim = apply_all_signatures(CSpGTCSm_0205_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pain_cog_emo');
                                           
appCSpGTCSm0205_emotion_sim = apply_all_signatures(CSpGTCSm_0205_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'kragelemotion');

appCSpGTCSm0205_pines_sim = apply_all_signatures(CSpGTCSm_0205_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pines');       
                                           
appCSpGTCSm0205_stroop_sim = apply_all_signatures(CSpGTCSm_0205_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'stroop');                                            
%% Passive Learning/Heterogeneous Sample, 0187tk
% similarity of CSpGTCSm contrast data and fear conditioning signature pattern
appCSpGTCSm0187_meta_fear_sim = apply_mask(CSpGTCSm_0187_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
CSplus_0187_meta_fear_expression = apply_mask(CSplus_0187,replace_empty(meta_fearcond_obj),  'pattern_expression', metric);
CSminus_0187_meta_fear_expression = apply_mask(CSminus_0187,replace_empty(meta_fearcond_obj),  'pattern_expression', metric);

% similarity of data masked with ROIs with pattern
nacc_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_nacc, replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
caudate_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_caudate, replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
amygdala_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_amygdala,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
cerebellum_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_cerebellum,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
nacc_CSplus_0187_meta_fear_expression = apply_mask(roi_CSplus_0187_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
nacc_CSminus_0187_meta_fear_expression = apply_mask(roi_CSminus_0187_nacc,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSplus_0187_meta_fear_expression = apply_mask(roi_CSplus_0187_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
caudate_CSminus_0187_meta_fear_expression = apply_mask(roi_CSminus_0187_caudate,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSplus_0187_meta_fear_expression = apply_mask(roi_CSplus_0187_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
amygdala_CSminus_0187_meta_fear_expression = apply_mask(roi_CSminus_0187_amygdala,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSplus_0187_meta_fear_expression = apply_mask(roi_CSplus_0187_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);
cerebellum_CSminus_0187_meta_fear_expression = apply_mask(roi_CSminus_0187_cerebellum,replace_empty(meta_fearcond_obj),  'pattern_expression', 'ignore_missing', metric);

% negative control patterns
appCSpGTCSm0187_cog_sim = apply_all_signatures(CSpGTCSm_0187_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pain_cog_emo');
                                           
appCSpGTCSm0187_emotion_sim = apply_all_signatures(CSpGTCSm_0187_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'kragelemotion');

appCSpGTCSm0187_pines_sim = apply_all_signatures(CSpGTCSm_0187_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'pines');                                       

appCSpGTCSm0187_stroop_sim = apply_all_signatures(CSpGTCSm_0187_con_obj, 'similarity_metric', metric,... 
                                               'conditionnames', {'appCSplusGTCSminus'},...
                                               'image_set', 'stroop'); 
                                           
%% -- visualize pattern and sample data --
%% whole brain fear conditioning pattern
% one subcortical cutaway
create_figure('fear conditoning pattern - whole brain')
my_display_obj = addbrain('right_cutaway');
my_display_obj = [my_display_obj addbrain('cerebellum')];
set(my_display_obj, 'FaceAlpha', 1);
render_on_surface(meta_fearcond_obj, my_display_obj);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% other cutaways
% possible keywords: {'left_cutaway' 'right_cutaway' 'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5'}
keywords = {'left_cutaway', 'accumbens_slab' 'coronal_slabs_4'};
for i = 1:length(keywords)
    create_figure('cutaways'); axis off
    surface_handles = surface(meta_fearcond_obj, keywords{i});
    drawnow, snapnow
end
clf()

%% contrast maps for each sample
whmontage = 5; 
plugin_check_or_create_slice_display;
allcons = {CSpGTCSm_0188_con_obj, CSpGTCSm_0205_con_obj, CSpGTCSm_0187_con_obj};     
k = length(allcons);
contrast_t_fdr = {};

o2 = removeblobs(o2);

for i = 1:k
    
    figtitle = sprintf('%s_05_scaleddata_FDR', allcons{i}.image_names{1});
    figstr = format_strings_for_legend(figtitle); 
    figstr = figstr{1};
    disp(figstr);

    contrast_t_fdr{i} = ttest(allcons{i}, .05, 'fdr');
    
    % 1st plot at 0.05 FDR
    % -----------------------------------------------
    o2 = removeblobs(o2);
    o2 = addblobs(o2, region(contrast_t_fdr{i}), 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]});
    
    axes(o2.montage{whmontage}.axis_handles(5));
    title(figstr, 'FontSize', 16)
    
    o2 = legend(o2);
    disp('colormap range of blobs in case legend gets messed up ')
    disp(o2.activation_maps{1}.cmaprange)
    drawnow, snapnow
    
end
%% -- Main Results -- 
%% whole brain
% Similarity of all sample contrasts to pattern
create_figure('Similarity of all sample contrasts to pattern - whole brain')
barplot_columns({appCSpGTCSm0188_meta_fear_sim
                 appCSpGTCSm0205_meta_fear_sim
                 appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% barplot all separate condition similarities to pattern
create_figure('Similarity of all sample conditions to pattern')
barplot_columns({CSplus_0188_meta_fear_expression 
                 CSminus_0188_meta_fear_expression 
                 CSplus_0205_meta_fear_expression 
                 CSminus_0205_meta_fear_expression 
                 CSplus_0187_meta_fear_expression
                 CSminus_0187_meta_fear_expression},...
                'names',{'appCS+' 
                         'appCS-'
                         'appCS+' 
                         'appCS-'
                         'appCS+' 
                         'appCS-'},...
                'colors',{color_CSpGTCSm_0188
                          color_CSpGTCSm_0188  
                          color_CSpGTCSm_0205
                          color_CSpGTCSm_0205
                          color_CSpGTCSm_0187
                          color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('conditions')
drawnow, snapnow

% ROC plot, classification of CSplus vs CSminus by pattern expression in all samples
create_figure('ROC plot all samples')
ROC_0188_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0188_meta_fear_expression; CSminus_0188_meta_fear_expression],...
                                               [ones(29,1);zeros(29,1)],...
                                               'twochoice','color',color_CSpGTCSm_0188,...
                                               'threshold_type', 'Optimal balanced error rate');
ROC_0205_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0205_meta_fear_expression; CSminus_0205_meta_fear_expression],...
                                               [ones(76,1);zeros(76,1)],...
                                                'twochoice','color',color_CSpGTCSm_0205,...
                                                'threshold_type', 'Optimal balanced error rate');
ROC_0187_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0187_meta_fear_expression; CSminus_0187_meta_fear_expression],...
                                               [ones(38,1);zeros(38,1)],...
                                                'twochoice','color',color_CSpGTCSm_0187,...
                                                'threshold_type', 'Optimal balanced error rate');
%ROC_legend = makelegend({'Active Hom', 'Active Het', 'Passive Het'}, {color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187});
drawnow, snapnow

%% NAcc
% Similarity of all sample contrasts to pattern - NAcc ROI 
create_figure('Similarity of all sample contrasts to pattern - NAcc')
barplot_columns({nacc_appCSpGTCSm0188_meta_fear_sim
                 nacc_appCSpGTCSm0205_meta_fear_sim
                 nacc_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% display pattern masked with NAcc
create_figure('fear conditoning pattern - NAcc'); axis off;
nacc_display_obj = addbrain('CIT168'); % cit168 contains nacc, caudate and put
nacc_display_obj = [nacc_display_obj addbrain('amygdala')];
nacc_display_obj = [nacc_display_obj addbrain('hippocampus')];
nacc_display_obj = [nacc_display_obj addbrain('brainstem')];
nacc_display_obj = [nacc_display_obj addbrain('cerebellum')];
set(nacc_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_nacc, nacc_display_obj);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

%% Caudate Nucleus
% Similarity of all sample contrasts to pattern - Caudate Nucleus ROI 
create_figure('Similarity of all sample contrasts to pattern - Caudate Nucleus')
barplot_columns({caudate_appCSpGTCSm0188_meta_fear_sim
                 caudate_appCSpGTCSm0205_meta_fear_sim
                 caudate_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% display pattern masked with Caudate Nucleus
create_figure('fear conditoning pattern - Caudate'); axis off;
caudate_display_obj = addbrain('CIT168');
caudate_display_obj = [caudate_display_obj addbrain('amygdala')];
caudate_display_obj = [caudate_display_obj addbrain('hippocampus')];
caudate_display_obj = [caudate_display_obj addbrain('brainstem')];
caudate_display_obj = [caudate_display_obj addbrain('cerebellum')];
set(caudate_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_caudate, caudate_display_obj);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

%% Amygdala
% Similarity of all sample contrasts to pattern - Amygdala ROI 
create_figure('Similarity of all sample contrasts to pattern - Amygdala')
barplot_columns({amygdala_appCSpGTCSm0188_meta_fear_sim
                 amygdala_appCSpGTCSm0205_meta_fear_sim
                 amygdala_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% display pattern masked with Amygdala
create_figure('fear conditoning pattern - Amygdala'); axis off;
amy_display_obj = addbrain('CIT168');
amy_display_obj = [amy_display_obj addbrain('amygdala')];
amy_display_obj = [amy_display_obj addbrain('hippocampus')];
amy_display_obj = [amy_display_obj addbrain('brainstem')];
amy_display_obj = [amy_display_obj addbrain('cerebellum')];
set(amy_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_amygdala, amy_display_obj);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

%% Cerebellum
% Similarity of all sample contrasts to pattern - Cerebellum ROI
create_figure('Similarity of all sample contrasts to pattern - Cerebellum')
barplot_columns({cerebellum_appCSpGTCSm0188_meta_fear_sim
                 cerebellum_appCSpGTCSm0205_meta_fear_sim
                 cerebellum_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% display pattern masked with Cerebellum
create_figure('fear conditoning pattern - Cerebellum'); axis off
cblm_display_obj = addbrain('cerebellum');
cblm_display_obj = [cblm_display_obj addbrain('CIT168')];
cblm_display_obj = [cblm_display_obj addbrain('amygdala')];
cblm_display_obj = [cblm_display_obj addbrain('hippocampus')];
cblm_display_obj = [cblm_display_obj addbrain('brainstem')];
set(cblm_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_cerebellum, cblm_display_obj);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow


%% -- Secondary Analysis results --
%% negative control - signatures
% check similarity with emotion signatures (Kragel 2015)
create_figure(strcat('Similarity of all sample contrasts to Fearful pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_emotion_sim.Fearful{:,:}
                 appCSpGTCSm0205_emotion_sim.Fearful{:,:}
                 appCSpGTCSm0187_emotion_sim.Fearful{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

create_figure(strcat('Similarity of all sample contrasts to Surprised pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_emotion_sim.Surprised{:,:}
                 appCSpGTCSm0205_emotion_sim.Surprised{:,:}
                 appCSpGTCSm0187_emotion_sim.Surprised{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

% Similarity of all sample contrasts to cog control pattern (Kragel 2018)
create_figure(strcat('Similarity of all sample contrasts to Cognitive Control pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_cog_sim.Cog_Wholebrain{:,:}
                 appCSpGTCSm0205_cog_sim.Cog_Wholebrain{:,:}
                 appCSpGTCSm0187_cog_sim.Cog_Wholebrain{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

% Similarity of all sample contrasts to negative affect pattern (Chang 2015)
create_figure(strcat('Similarity of all sample contrasts to PINES pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_pines_sim.PINES{:,:}
                 appCSpGTCSm0205_pines_sim.PINES{:,:}
                 appCSpGTCSm0187_pines_sim.PINES{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

% Similarity of all sample contrasts to Stroop pattern (Silvestrini 2020)
create_figure(strcat('Similarity of all sample contrasts to stroop pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_stroop_sim.Stroop{:,:}
                 appCSpGTCSm0205_stroop_sim.Stroop{:,:}
                 appCSpGTCSm0187_stroop_sim.Stroop{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow


%% fear conditioning pattern - table of regions
% Table of results unc 
r_meta_fearcond = region(meta_fearcond_obj);    
[rmetafearpos, rmetafearneg] = table(r_meta_fearcond);       % add labels
drawnow, snapnow

%% SVM appetitive CSplusGTCSminus predictive weight maps
%% ROC plot
% [dist_from_hyperplane_0188, Y_0188, svm_dist_pos_neg_0188, svm_dist_pos_neg_matrix_0188] = plugin_svm_contrasts_get_results_per_subject(MID_0188_DAT, MID_0188_svm_stats_results, MID_0188_DATA_OBJ);
% [dist_from_hyperplane_0205, Y_0205, svm_dist_pos_neg_0205, svm_dist_pos_neg_matrix_0205] = plugin_svm_contrasts_get_results_per_subject(MID_0205_DAT, MID_0205_svm_stats_results, MID_0205_DATA_OBJ);
% [dist_from_hyperplane_0187, Y_0187, svm_dist_pos_neg_0187, svm_dist_pos_neg_matrix_0187] = plugin_svm_contrasts_get_results_per_subject(Acq_0187_DAT, Acq_0187_svm_stats_results, Acq_0187_DATA_OBJ);
% 
% create_figure('SVM appCS+ > appCS- classifier - ROC plot all samples')
% ROC_0188_CSp_CSm_SVM_forced_choice = roc_plot(dist_from_hyperplane_0188{1}, logical(Y_0188{1} > 0), 'color', color_CSpGTCSm_0188, 'twochoice');
% ROC_0205_CSp_CSm_SVM_forced_choice = roc_plot(dist_from_hyperplane_0205{1}, logical(Y_0205{1} > 0), 'color', color_CSpGTCSm_0205, 'twochoice');
% ROC_0187_CSp_CSm_SVM_forced_choice = roc_plot(dist_from_hyperplane_0187{1}, logical(Y_0187{1} > 0), 'color', color_CSpGTCSm_0187, 'twochoice');
% ROC_legend = makelegend({'Active Hom', 'Active Het', 'Passive Het'}, {color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187});
% drawnow, snapnow

%% tables of region for all samples 
CSpGTCSm_0188_svm_obj_fdr_05 = threshold(CSpGTCSm_0188_svm_obj, 0.05, 'fdr');
CSpGTCSm_0205_svm_obj_fdr_05 = threshold(CSpGTCSm_0205_svm_obj, 0.05, 'fdr');
CSpGTCSm_0187_svm_obj_fdr_05 = threshold(CSpGTCSm_0187_svm_obj, 0.05, 'fdr');

%% Active Learning/Homogeneous Sample regions
% fdr .05 corr
r_CSpGTCSm_0188_svm_fdr_05 = region(CSpGTCSm_0188_svm_obj_fdr_05);    
    [rCSp0188_svm_fdr_05, rCSm0188_svm_fdr_05] = table(r_CSpGTCSm_0188_svm_fdr_05);       % add labels
drawnow, snapnow

% unc .05
r_CSpGTCSm_0188_svm = region(CSpGTCSm_0188_svm_obj);    
    [rCSp0188_svm, rCSm0188_svm] = table(r_CSpGTCSm_0188_svm);       % add labels
drawnow, snapnow

%T = table(r_CSpGTCSm_0188_svm);
%writetable(T, fullfile(base_dir, '{scripts/html/r_CSpGTCSm_0188_svm.xlsx'))
%% Active Learning/Heterogeneous Sample regions
% fdr .05 corr
r_CSpGTCSm_0205_svm_fdr_05 = region(CSpGTCSm_0205_svm_obj_fdr_05);    
    [rCSp0205_svm_fdr_05, rCSm0205_svm_fdr_05] = table(r_CSpGTCSm_0205_svm_fdr_05);       % add labels
drawnow, snapnow

% unc .05
r_CSpGTCSm_0205_svm = region(CSpGTCSm_0205_svm_obj);    
    [rCSp0205_svm, rCSm0205_svm] = table(r_CSpGTCSm_0205_svm);       % add labels
drawnow, snapnow

%% Passive Learning/Heterogeneous Sample regions
% fdr .05 corr
r_CSpGTCSm_0187_svm_fdr_05 = region(CSpGTCSm_0187_svm_obj_fdr_05);    
    [rCSp0187_svm_fdr_05, rCSm0187_svm_fdr_05] = table(r_CSpGTCSm_0187_svm_fdr_05);       % add labels
drawnow, snapnow

% unc .05
r_CSpGTCSm_0187_svm = region(CSpGTCSm_0187_svm_obj);    
    [rCSp0187_svm, rCSm0187_svm] = table(r_CSpGTCSm_0187_svm);       % add labels
drawnow, snapnow


%% display appetitive CSplusGTCSminus predictive weight maps 
%% sample 1
display_0188_svm_2 = orthviews(CSpGTCSm_0188_svm_obj);
drawnow, snapnow
clf();
display_obj_0188_svm = canlab_results_fmridisplay(CSpGTCSm_0188_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow

%% sample 2
display_0205_svm_2 = orthviews(CSpGTCSm_0205_svm_obj);
drawnow, snapnow
clf();
display_obj_0205_svm = canlab_results_fmridisplay(CSpGTCSm_0205_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow

%% sample 3
display_0187_svm_2 = orthviews(CSpGTCSm_0187_svm_obj);
drawnow, snapnow
clf();
display_obj_0187_svm = canlab_results_fmridisplay(CSpGTCSm_0187_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow


