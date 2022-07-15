%% author: s-kline
% analysis for paper comparing appetitive and aversive conditioning
clear all, close all
%% -- prep -- 
%% set up directories
base_dir = '\aversive appetitive conditioning\github_repo\';
% multivariate pattern dir 
meta_dir = fullfile(base_dir, 'mv patterns\n_signature_fear_conditioning\');

% appetitive cond sample dirs
MID_0188_dir = fullfile(base_dir,'cond samples\0188ok\CSplus_classifier_onlyPlacebo_for_paper\'); % active hom sample
MID_0205_dir = fullfile(base_dir,'cond samples\0205ok\0205_t1_CSplus_classifier_for_paper\'); % active het sample
Acq_0187_dir = fullfile(base_dir,'cond samples\0187tk\0187_CSplus_classifier_onlyaware_for_paper\'); % passive het sample

%% load data
% -- neural fear conditiong signature map 
% meta analysis, Fullana et al 2016 
meta_fearcond_obj = statistic_image(fullfile(meta_dir, 'MyMean_z.nii.gz')); % map of activations and deactivations related to aversive CS+
color_meta_fearcond = [224  0  97]/255; % ruby (E00061)

% create sig-masked image with both positive and negative values
meta_fearcond_pos_05 = statistic_image(fullfile(meta_dir, 'Fullana_2016_Figure1.nii.gz'));
meta_fearcond_neg_05 = statistic_image(fullfile(meta_dir, 'Fullana_2016_Figure2.nii.gz'));
mpos_05 = fmri_mask_image(meta_fearcond_pos_05);
mneg_05 = fmri_mask_image(meta_fearcond_neg_05);
meta_fear_thresh_05_mask = image_math(mpos_05, mneg_05, 'add');
meta_fearcond_thr05_obj = apply_mask(meta_fearcond_obj, meta_fear_thresh_05_mask);

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
% Putamen
putamen = select_atlas_subset(basal_ganglia, 'labels', {'Putamen_Pa_L', 'Putamen_Pa_R', 'Putamen_Pp_L', 'Putamen_Pp_R'}, 'flatten');
putamen = resample_space(putamen, CSpGTCSm_0188_con_obj);
putamen.fullpath='putamen.nii';putamen.write('overwrite');
% Thalamus
thalamus = fmri_mask_image('Thalamus_thr50_bilateral.nii'); % from Harvard Oxford Cortical, no need to resample
% Insula
insula = fmri_mask_image('Insula_thr50_bilateral.nii'); % from Harvard Oxford Cortical, no need to resample

% -------------------------------------------------------------------------


%% mask sample data
% Active Learning/Homogeneous Sample 
roi_CSpGTCSm_0188_con_obj_nacc = apply_mask(CSpGTCSm_0188_con_obj, nacc);
roi_CSpGTCSm_0188_con_obj_caudate = apply_mask(CSpGTCSm_0188_con_obj, caudate);
roi_CSpGTCSm_0188_con_obj_amygdala = apply_mask(CSpGTCSm_0188_con_obj, amygdala);
roi_CSpGTCSm_0188_con_obj_cerebellum = apply_mask(CSpGTCSm_0188_con_obj, cerebellum);
roi_CSpGTCSm_0188_con_obj_putamen = apply_mask(CSpGTCSm_0188_con_obj, putamen);
roi_CSpGTCSm_0188_con_obj_thalamus = apply_mask(CSpGTCSm_0188_con_obj, thalamus);
roi_CSpGTCSm_0188_con_obj_insula = apply_mask(CSpGTCSm_0188_con_obj, insula);
% -------------------------------------------------------------------------
% Active Learning/Heterogeneous Sample
roi_CSpGTCSm_0205_con_obj_nacc = apply_mask(CSpGTCSm_0205_con_obj, nacc);
roi_CSpGTCSm_0205_con_obj_caudate = apply_mask(CSpGTCSm_0205_con_obj, caudate);
roi_CSpGTCSm_0205_con_obj_amygdala = apply_mask(CSpGTCSm_0205_con_obj, amygdala);
roi_CSpGTCSm_0205_con_obj_cerebellum = apply_mask(CSpGTCSm_0205_con_obj, cerebellum);
roi_CSpGTCSm_0205_con_obj_putamen = apply_mask(CSpGTCSm_0205_con_obj, putamen);
roi_CSpGTCSm_0205_con_obj_thalamus = apply_mask(CSpGTCSm_0205_con_obj, thalamus);
roi_CSpGTCSm_0205_con_obj_insula = apply_mask(CSpGTCSm_0205_con_obj, insula);
% -------------------------------------------------------------------------
% Passive Learning/Heterogeneous Sample
roi_CSpGTCSm_0187_con_obj_nacc = apply_mask(CSpGTCSm_0187_con_obj, nacc);
roi_CSpGTCSm_0187_con_obj_caudate = apply_mask(CSpGTCSm_0187_con_obj, caudate);
roi_CSpGTCSm_0187_con_obj_amygdala = apply_mask(CSpGTCSm_0187_con_obj, amygdala);
roi_CSpGTCSm_0187_con_obj_cerebellum = apply_mask(CSpGTCSm_0187_con_obj, cerebellum);
roi_CSpGTCSm_0187_con_obj_putamen = apply_mask(CSpGTCSm_0187_con_obj, putamen);
roi_CSpGTCSm_0187_con_obj_thalamus = apply_mask(CSpGTCSm_0187_con_obj, thalamus);
roi_CSpGTCSm_0187_con_obj_insula = apply_mask(CSpGTCSm_0187_con_obj, insula);


% -------------------------------------------------------------------------
%% -- Main Analysis -- 
%% test meta fear conditioning signature on appetitive cond data -- 
% ROI analysis: masked data with ROIs and calculated similarity
% the canlab_pattern_similarity doc implies, that the pattern should not be
% masked, because of how the function deals with zeroes in statistic images vs data
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
putamen_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_putamen,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
thalamus_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_thalamus,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
insula_appCSpGTCSm0188_meta_fear_sim = apply_mask(roi_CSpGTCSm_0188_con_obj_insula,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);

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
putamen_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_putamen,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
thalamus_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_thalamus,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
insula_appCSpGTCSm0205_meta_fear_sim = apply_mask(roi_CSpGTCSm_0205_con_obj_insula,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);

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
putamen_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_putamen,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
thalamus_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_thalamus,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);
insula_appCSpGTCSm0187_meta_fear_sim = apply_mask(roi_CSpGTCSm_0187_con_obj_insula,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric);


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
% define some stuff for all figures 
% same color values for all activation maps
min_neg_col = [0 0 1];
max_neg_col = [0 1 1];
min_pos_col = [1 .5 0];
max_pos_col = [1 1 0];
pos_colormap = colormap_tor(min_pos_col, max_pos_col);
neg_colormap = colormap_tor(min_neg_col, max_neg_col);
% for lines on the saggital slices
linewidth = 4
% values where the map was thresholded
vals = unique(meta_fearcond_thr05_obj.dat);
neg_lim = max(vals(vals<0)); % largest negative value
pos_lim = min(vals(vals>0)); % smallest positive value


% cutaways used in fig. 2
keywords = {'left_cutaway', 'accumbens_slab' 'coronal_slabs_4'};
%keywords = {'coronal_slabs_4'};
for i = 1:length(keywords)
    create_figure('cutaways'); axis off
    surface_handles = surface(meta_fearcond_thr05_obj, keywords{i}, 'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap);
    drawnow, snapnow
end
clf()

% one subcortical cutaway, used in fig. 3
create_figure('fear conditoning pattern - whole brain')
my_display_obj = addbrain('right_cutaway');
my_display_obj = [my_display_obj addbrain('cerebellum')];
set(my_display_obj, 'FaceAlpha', 1); % make surface solid
render_on_surface(meta_fearcond_thr05_obj, my_display_obj,... 
                'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% wholebrain slice display used in fig 3
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2,  meta_fearcond_thr05_obj,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-3.65 9.7])
    %'cmaprange', [-3.65 -1.31 3.36 9.7]
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = 0;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', 3);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh1);
% o2 = removeblobs(slices);
o2 = addblobs(o2,  meta_fearcond_thr05_obj,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-3.65 9.7])
    %'cmaprange', [-3.65 -1.31 3.36 9.7]);
% make a legend with a white box to indicate where the map was thresholded
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
drawnow, snapnow

% mask pattern (only for roi figures)
roi_meta_fearcond_obj_nacc = apply_mask(meta_fearcond_thr05_obj, nacc);
roi_meta_fearcond_obj_caudate = apply_mask(meta_fearcond_thr05_obj, caudate);
roi_meta_fearcond_obj_amygdala = apply_mask(meta_fearcond_thr05_obj, amygdala);
roi_meta_fearcond_obj_cerebellum = apply_mask(meta_fearcond_thr05_obj, cerebellum);
roi_meta_fearcond_obj_putamen = apply_mask(meta_fearcond_thr05_obj, putamen);
roi_meta_fearcond_obj_thalamus = apply_mask(meta_fearcond_thr05_obj, thalamus);
roi_meta_fearcond_obj_insula = apply_mask(meta_fearcond_thr05_obj, insula);

%% contrast maps for each sample
whmontage = 5; 
%plugin_check_or_create_slice_display;
create_figure('fmridisplay'); axis off
o2 = canlab_results_fmridisplay([], 'noverbose', 'compact');
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
    
    % plot at 0.05 FDR
    % -----------------------------------------------
    o2 = removeblobs(o2);
    o2 = addblobs(o2, region(contrast_t_fdr{i}),'trans',...
        'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col});
    
    o2 = addblobs(o2, nacc, 'color', [0.9 0.9 1]); % black
    %o2 = addblobs(o2, caudate, 'outline', 'color', [0 0.5 0]); % dark green
    
    axes(o2.montage{whmontage}.axis_handles(5));
    title(figstr, 'FontSize', 16)
    
    o2 = legend(o2);
    disp('colormap range of blobs in case legend gets messed up ')
    disp(o2.activation_maps{1}.cmaprange)
    drawnow, snapnow
    
end

% draw another montage with conditioning regions as blobs
title('anatomical regions', 'FontSize', 16)
o2 = removeblobs(o2);
o2 = addblobs(o2, nacc, 'outline', 'color', [0 1 0]);
drawnow, snapnow

o2 = removeblobs(o2);
o2 = addblobs(o2, caudate, 'outline', 'color', [0 0.5 1]);
drawnow, snapnow


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

% fig 3: pattern masked with NAcc, 3D surface
create_figure('fear conditoning pattern - NAcc'); axis off;
nacc_display_obj = addbrain('CIT168'); % cit168 contains nacc, caudate and put
nacc_display_obj = [nacc_display_obj addbrain('amygdala')];
nacc_display_obj = [nacc_display_obj addbrain('hippocampus')];
nacc_display_obj = [nacc_display_obj addbrain('brainstem')];
nacc_display_obj = [nacc_display_obj addbrain('cerebellum')];
%nacc_display_obj = [nacc_display_obj addbrain('hires right')];
set(nacc_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_nacc, nacc_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-0.46 6.26]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with NAcc, ax and sag slices
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_nacc,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-0.46 6.26]...
    );
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = -10;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
% o2 = removeblobs(slices);
o2 = addblobs(o2,  roi_meta_fearcond_obj_nacc,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-0.46 6.26]...
    );
% make a legend with a white box to indicate where the map was thresholded
% make a legend with a white box to indicate where the map was thresholded
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
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

% fig 3: pattern masked with Caudate Nucleus, 3D surface
create_figure('fear conditoning pattern - Caudate'); axis off;
caudate_display_obj = addbrain('CIT168');
caudate_display_obj = [caudate_display_obj addbrain('amygdala')];
caudate_display_obj = [caudate_display_obj addbrain('hippocampus')];
caudate_display_obj = [caudate_display_obj addbrain('brainstem')];
caudate_display_obj = [caudate_display_obj addbrain('cerebellum')];
%caudate_display_obj = [caudate_display_obj addbrain('hires right')];
set(caudate_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_caudate, caudate_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-1.83 7.82]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with caudate, ax slice
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_caudate,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.83 7.82]...
    );
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = 2;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
% o2 = removeblobs(slices);
o2 = addblobs(o2,  roi_meta_fearcond_obj_caudate,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.83 7.82]...
    );
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
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

% fig 3: pattern masked with Amygdala, 3D surface
create_figure('fear conditoning pattern - Amygdala'); axis off;
amy_display_obj = addbrain('CIT168');
amy_display_obj = [amy_display_obj addbrain('amygdala')];
amy_display_obj = [amy_display_obj addbrain('hippocampus')];
amy_display_obj = [amy_display_obj addbrain('brainstem')];
amy_display_obj = [amy_display_obj addbrain('cerebellum')];
%amy_display_obj = [amy_display_obj addbrain('hires right')];
set(amy_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_amygdala, amy_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-1.61 5.46]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with amy, ax slice
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_amygdala,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.61 5.46]...
    );
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = -16;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
% o2 = removeblobs(slices);
o2 = addblobs(o2,  roi_meta_fearcond_obj_amygdala,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.61 5.46]...
    );
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
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

% fig 3: pattern masked with cerebellum, 3D surface
create_figure('fear conditoning pattern - Cerebellum'); axis off
cblm_display_obj = addbrain('cerebellum');
cblm_display_obj = [cblm_display_obj addbrain('CIT168')];
cblm_display_obj = [cblm_display_obj addbrain('amygdala')];
cblm_display_obj = [cblm_display_obj addbrain('hippocampus')];
cblm_display_obj = [cblm_display_obj addbrain('brainstem')];
%cblm_display_obj = [cblm_display_obj addbrain('hires right')];
set(cblm_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_cerebellum, cblm_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-2.51 5.14]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow
clf()

% fig 3: pattern masked with cerebellum, ax and sag slices
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_cerebellum,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.51 5.14]...
    );
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = -28;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
% o2 = removeblobs(slices);
o2 = addblobs(o2,  roi_meta_fearcond_obj_cerebellum,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.51 5.14]...
    );
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
drawnow, snapnow
%% Putamen
% Similarity of all sample contrasts to pattern - Putamen ROI
create_figure('Similarity of all sample contrasts to pattern - Putamen')
barplot_columns({putamen_appCSpGTCSm0188_meta_fear_sim
                 putamen_appCSpGTCSm0205_meta_fear_sim
                 putamen_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% fig 3: pattern masked with putamen, 3D surface
create_figure('fear conditoning pattern - Putamen'); axis off
put_display_obj = addbrain('CIT168');
put_display_obj = [put_display_obj addbrain('amygdala')];
put_display_obj = [put_display_obj addbrain('hippocampus')];
put_display_obj = [put_display_obj addbrain('brainstem')];
put_display_obj = [put_display_obj addbrain('cerebellum')];
%amy_display_obj = [amy_display_obj addbrain('hires right')];
set(put_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_putamen, put_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-2.56 7.97]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with putamen, ax and sag slices
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_putamen,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.56 7.97]);
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = 2;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
o2 = addblobs(o2,  roi_meta_fearcond_obj_putamen,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.56 7.97]);
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
drawnow, snapnow

%% Thalamus
% Similarity of all sample contrasts to pattern - Thalamus ROI
create_figure('Similarity of all sample contrasts to pattern - Thalamus')
barplot_columns({thalamus_appCSpGTCSm0188_meta_fear_sim
                 thalamus_appCSpGTCSm0205_meta_fear_sim
                 thalamus_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% fig 3: pattern masked with thalamus, 3D surface
create_figure('fear conditoning pattern - Thalamus'); axis off
thal_display_obj = addbrain('CIT168');
thal_display_obj = [thal_display_obj addbrain('amygdala')];
thal_display_obj = [thal_display_obj addbrain('hippocampus')];
thal_display_obj = [thal_display_obj addbrain('brainstem')];
thal_display_obj = [thal_display_obj addbrain('cerebellum')];
set(thal_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_thalamus, thal_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-1.49 6.65]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with thalamus, ax and sag slices
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_thalamus,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.49 6.65]);
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = 7;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
o2 = addblobs(o2,  roi_meta_fearcond_obj_thalamus,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-1.49 6.65]);
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
drawnow, snapnow

%% Insula
% Similarity of all sample contrasts to pattern - Insula ROI
create_figure('Similarity of all sample contrasts to pattern - Insula')
barplot_columns({insula_appCSpGTCSm0188_meta_fear_sim
                 insula_appCSpGTCSm0205_meta_fear_sim
                 insula_appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
ylim([-1.5 1.5])
drawnow, snapnow

% fig 3: pattern masked with insula, 3D surface
create_figure('fear conditoning pattern - Insula'); axis off
ins_display_obj = addbrain('CIT168');
ins_display_obj = [ins_display_obj addbrain('amygdala')];
ins_display_obj = [ins_display_obj addbrain('hippocampus')];
ins_display_obj = [ins_display_obj addbrain('brainstem')];
ins_display_obj = [ins_display_obj addbrain('cerebellum')];
set(ins_display_obj, 'FaceAlpha', 1);
render_on_surface(roi_meta_fearcond_obj_insula, ins_display_obj,...
    'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap,...
    'clim', [-2.58 9.7]);
view(222, 5);
lightRestoreSingle;
drawnow, snapnow

% fig 3: pattern masked with insula, ax and sag slices
clf()
slices = fmridisplay();
% saggital slice
axh2 = axes('Position', [-0.1 0.2 .7 .7]);
o2 = montage(slices, 'sagittal', 'wh_slice', [0 0 0], 'onerow', 'existing_axes', axh2);
o2 = removeblobs(o2);
o2 = addblobs(o2, roi_meta_fearcond_obj_insula,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.58 9.7]);
% Line to indicate axial slice
axes(o2.montage{1}.axis_handles)
loc = 0;
hh(1) = plot([-125 85], [loc loc], 'b', 'LineWidth', linewidth);
% axial slice
axh1 = axes('Position', [0.3 0.2 .7 .7]);
o2 = montage(slices, 'axial', 'wh_slice', [0 0 loc], 'onerow', 'existing_axes', axh1);
o2 = addblobs(o2,  roi_meta_fearcond_obj_insula,...
    'splitcolor', {min_neg_col max_neg_col min_pos_col max_pos_col},...
    'cmaprange', [-2.58 9.7]);
o2 = legend(o2, 'newfig')
plot([neg_lim pos_lim], [0.5 0.5], 'white', 'Linewidth', 13)
plot([0 0], [0 1], 'black', 'Linewidth', 1) % indicate zero
drawnow, snapnow

%% -- Secondary Analysis results --
%% negative control - emotion signatures
% check similarity with emotion signatures
for i = 1:length(appCSpGTCSm0188_emotion_sim.signaturenames)
    sig = char(appCSpGTCSm0188_emotion_sim.signaturenames(i));
  
    sig_sim_0188 = getfield(appCSpGTCSm0188_emotion_sim, sig);
    sig_sim_0205 = getfield(appCSpGTCSm0205_emotion_sim, sig);
    sig_sim_0187 = getfield(appCSpGTCSm0187_emotion_sim, sig);
    
    % Similarity of all sample contrasts to current signature
    create_figure(strcat('Similarity of all sample contrasts to ', sig, ' - whole brain'))
    barplot_columns({sig_sim_0188{:,:}
                     sig_sim_0205{:,:}
                     sig_sim_0187{:,:}},...
                    'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                    'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
    ylabel('Cosine Similarity')
    xlabel('appCS+ > appCS- contrasts')
    %ylim([-1.5 1.5])
    drawnow, snapnow
    
    % test if meta fear pattern has higher mean sim values than control:
    [h,p,ci,stats] = ttest(appCSpGTCSm0188_meta_fear_sim, sig_sim_0188{:,:})
    [h,p,ci,stats] = ttest(appCSpGTCSm0205_meta_fear_sim, sig_sim_0205{:,:})
    [h,p,ci,stats] = ttest(appCSpGTCSm0187_meta_fear_sim, sig_sim_0187{:,:})
end

%% negative control - cognitive control pattern
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

% test if meta fear pattern has higher mean sim values than control:
[h,p,ci,stats] = ttest(appCSpGTCSm0188_meta_fear_sim, appCSpGTCSm0188_cog_sim.Cog_Wholebrain{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0205_meta_fear_sim, appCSpGTCSm0205_cog_sim.Cog_Wholebrain{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0187_meta_fear_sim, appCSpGTCSm0187_cog_sim.Cog_Wholebrain{:,:})

%% negative control - negative affect pattern
% Similarity of all sample contrasts to negative affect pattern ()
create_figure(strcat('Similarity of all sample contrasts to PINES pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_pines_sim.PINES{:,:}
                 appCSpGTCSm0205_pines_sim.PINES{:,:}
                 appCSpGTCSm0187_pines_sim.PINES{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

% test if meta fear pattern has higher mean sim values than control:
[h,p,ci,stats] = ttest(appCSpGTCSm0188_meta_fear_sim, appCSpGTCSm0188_pines_sim.PINES{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0205_meta_fear_sim, appCSpGTCSm0205_pines_sim.PINES{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0187_meta_fear_sim, appCSpGTCSm0187_pines_sim.PINES{:,:})

%% negative control - stroop pattern
% Similarity of all sample contrasts to Stroop pattern ()
create_figure(strcat('Similarity of all sample contrasts to stroop pattern - whole brain'))
barplot_columns({appCSpGTCSm0188_stroop_sim.Stroop{:,:}
                 appCSpGTCSm0205_stroop_sim.Stroop{:,:}
                 appCSpGTCSm0187_stroop_sim.Stroop{:,:}},...
                'names',{'Active Hom', 'Active Het', 'Passive Het'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('appCS+ > appCS- contrasts')
drawnow, snapnow

[h,p,ci,stats] = ttest(appCSpGTCSm0188_meta_fear_sim, appCSpGTCSm0188_stroop_sim.Stroop{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0205_meta_fear_sim, appCSpGTCSm0205_stroop_sim.Stroop{:,:})
[h,p,ci,stats] = ttest(appCSpGTCSm0187_meta_fear_sim, appCSpGTCSm0187_stroop_sim.Stroop{:,:})

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

SVM_0188_meta_fear_sim = apply_mask(CSpGTCSm_0188_svm_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric)
SVM_0205_meta_fear_sim = apply_mask(CSpGTCSm_0205_svm_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric)
SVM_0187_meta_fear_sim = apply_mask(CSpGTCSm_0187_svm_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', metric)

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

