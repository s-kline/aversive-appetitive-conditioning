%% author: s-kline
% analysis for paper comparing appetitive and aversive conditioning
clear all, close all

%% compare aversive and appetitive conditioning
base_dir = '\aversive-appetitive-conditioning\';
% multivariate pattern dirs 
meta_dir = fullfile(base_dir, 'n_signature_fear_conditioning\');
% appetitive cond sample dirs
MID_0188_dir = fullfile(base_dir,'cond samples\ActiveLearning_Homogeneous\');
MID_0205_dir = fullfile(base_dir,'cond samples\ActiveLearning_Heterogeneous\');
Acq_0187_dir = fullfile(base_dir,'cond samples\PassiveLearning_Heterogeneous\');

%% load patterns of interest
% neural fear conditiong signature map (from meta analysis, Fullana et al 2015) 
meta_fearcond_obj = statistic_image(fullfile(meta_dir, 'MyMean_z.nii.gz')); % map of activations and deactivations related to aversive CS+
color_meta_fearcond = [224  0  97]/255; % ruby (E00061)

%% load appetitive conditioning data (0188ok study, monetary incentive delay) 
color_CSpGTCSm_0188 = [76 201 240]/255; % vivid sky blue (4CC9F06)

MID_0188_DATA_OBJ = importdata(fullfile(MID_0188_dir, '\results\data_objects.mat'));
MID_0188_CONTRASTS_OBJ = importdata(fullfile(MID_0188_dir, '\results\contrast_data_objects.mat'));
MID_0188_svm_stats_results = importdata(fullfile(MID_0188_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
CSplus_0188 = fmri_data(MID_0188_DATA_OBJ{1, 3});
CSminus_0188 = fmri_data(MID_0188_DATA_OBJ{1, 6});
CSpGTCSm_0188_con_obj = MID_0188_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1}; 

CSpGTCSm_0188_svm_obj = MID_0188_svm_stats_results{1,1}.weight_obj;
fname = 'CSpGTCSm_0188_svm_obj.nii';
write(CSpGTCSm_0188_svm_obj, 'fname', fullfile(base_dir, 'scripts\', fname), 'overwrite');

%% load appetitive conditioning data (0205ok study, monetary incentive delay)
color_CSpGTCSm_0205 = [58 12 163]/255; % blue (trypan blue 3A0CA3)

MID_0205_DATA_OBJ = importdata(fullfile(MID_0205_dir, '\results\data_objects.mat'));
MID_0205_CONTRASTS_OBJ = importdata(fullfile(MID_0205_dir, '\results\contrast_data_objects.mat'));
MID_0205_svm_stats_results = importdata(fullfile(MID_0205_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
CSplus_0205 = fmri_data(MID_0205_DATA_OBJ{1, 3});
CSminus_0205 = fmri_data(MID_0205_DATA_OBJ{1, 6});
CSpGTCSm_0205_con_obj = MID_0205_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1};

CSpGTCSm_0205_svm_obj = MID_0205_svm_stats_results{1,1}.weight_obj;
fname = 'CSpGTCSm_0205_svm_obj.nii';
write(CSpGTCSm_0205_svm_obj, 'fname', fullfile(base_dir, 'scripts\', fname), 'overwrite');

%% load appetitive conditioning data (0187tk study, passive taskk)
color_CSpGTCSm_0187 = [247 37 133]/255; % flickr pink (F72585)

Acq_0187_DATA_OBJ = importdata(fullfile(Acq_0187_dir, '\results\data_objects.mat'));
Acq_0187_CONTRASTS_OBJ = importdata(fullfile(Acq_0187_dir, '\results\contrast_data_objects.mat'));
Acq_0187_svm_stats_results = importdata(fullfile(Acq_0187_dir, '\results\svm_stats_results_contrasts_masked.mat'));
% load data objects and svm results
CSplus_0187 = fmri_data(Acq_0187_DATA_OBJ{1, 1});
CSminus_0187 = fmri_data(Acq_0187_DATA_OBJ{1, 2});
CSpGTCSm_0187_con_obj = Acq_0187_CONTRASTS_OBJ.DATA_OBJ_CONsc{1,1}; 

CSpGTCSm_0187_svm_obj = Acq_0187_svm_stats_results{1,1}.weight_obj;
fname = 'CSpGTCSm_0187_svm_obj.nii';
write(CSpGTCSm_0187_svm_obj, 'fname', fullfile(base_dir, 'scripts\', fname), 'overwrite');

%% Sample 1 -- 0188tk -- 
%% test meta fear conditioning signature on appetitive cond data (0188ok)
% similarity of CSpGTCSm classifier map and fear conditioning signature pattern
appCSpGTCSm0188_meta_fear_sim = apply_mask(CSpGTCSm_0188_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', 'cosine_similarity');
CSplus_0188ok_meta_fear_expression = apply_mask(CSplus_0188,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');
CSminus_0188ok_meta_fear_expression = apply_mask(CSminus_0188,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');

%% Sample 2 - 0205ok -- 
%% test meta fear conditioning signature on appetitive cond data (0205ok)
% similarity of CSpGTCSm classifier map and threat signature pattern
appCSpGTCSm0205_meta_fear_sim = apply_mask(CSpGTCSm_0205_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', 'cosine_similarity');
CSplus_0205ok_meta_fear_expression = apply_mask(CSplus_0205,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');
CSminus_0205ok_meta_fear_expression = apply_mask(CSminus_0205,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');

%% Sample 3 - 0187tk -- 
%% test meta fear conditioning signature on appetitive cond data (0187tk)
% similarity of CSpGTCSm classifier map and fear conditioning signature pattern
appCSpGTCSm0187_meta_fear_sim = apply_mask(CSpGTCSm_0187_con_obj,replace_empty(meta_fearcond_obj), 'pattern_expression', 'ignore_missing', 'cosine_similarity');
CSplus_0187tk_meta_fear_expression = apply_mask(CSplus_0187,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');
CSminus_0187tk_meta_fear_expression = apply_mask(CSminus_0187,replace_empty(meta_fearcond_obj),  'pattern_expression', 'cosine_similarity');

%% plots --
%% barplot all separate condition similarities to pattern
create_figure('Similarity of all sample conditions to pattern')
barplot_columns({CSplus_0188ok_meta_fear_expression 
                 CSminus_0188ok_meta_fear_expression 
                 CSplus_0205ok_meta_fear_expression 
                 CSminus_0205ok_meta_fear_expression 
                 CSplus_0187tk_meta_fear_expression
                 CSminus_0187tk_meta_fear_expression},...
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

%% barplot for all contrast similarities
create_figure('Similarity of all sample contrasts to pattern')
barplot_columns({appCSpGTCSm0188_meta_fear_sim
                 appCSpGTCSm0205_meta_fear_sim
                 appCSpGTCSm0187_meta_fear_sim},...
                'names',{'Sample 1', 'Sample 2', 'Sample 3'},...
                'colors',{color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187},'nofigure')
ylabel('Cosine Similarity')
xlabel('CS+ > CS- contrasts')
drawnow, snapnow

%% ROC plots for all samples 
create_figure('ROC plot all samples')
ROC_0188_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0188ok_meta_fear_expression; CSminus_0188ok_meta_fear_expression],...
                                               [ones(29,1);zeros(29,1)],...
                                               'twochoice','color',color_CSpGTCSm_0188,...
                                               'threshold_type', 'Optimal balanced error rate');
ROC_0205_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0205ok_meta_fear_expression; CSminus_0205ok_meta_fear_expression],...
                                               [ones(76,1);zeros(76,1)],...
                                                'twochoice','color',color_CSpGTCSm_0205,...
                                                'threshold_type', 'Optimal balanced error rate');
ROC_0187_CSp_CSm_meta_forced_choice = roc_plot([CSplus_0187tk_meta_fear_expression; CSminus_0187tk_meta_fear_expression],...
                                               [ones(38,1);zeros(38,1)],...
                                                'twochoice','color',color_CSpGTCSm_0187,...
                                                'threshold_type', 'Optimal balanced error rate');
ROC_legend = makelegend({'Sample 1', 'Sample 2', 'Sample 3'}, {color_CSpGTCSm_0188, color_CSpGTCSm_0205, color_CSpGTCSm_0187});
drawnow, snapnow

%% display MA fear conditioning predictive weight map with table of regions
% Table of results unc 
r_meta_fearcond = region(meta_fearcond_obj);    
[rmetafearpos, rmetafearneg] = table(r_meta_fearcond);       % add labels
drawnow, snapnow

% one subcortical cutaway
create_figure('MA pattern')
my_display_obj = surface(meta_fearcond_obj);
drawnow, snapnow

% other cutaways
% possible keywords: {'left_cutaway' 'right_cutaway' 'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5'}
keywords = {'accumbens_slab' 'coronal_slabs_4'};
for i = 1:length(keywords)
    create_figure('cutaways'); axis off
    surface_handles = surface(meta_fearcond_obj, keywords{i});
    drawnow, snapnow
end
clf()

display_obj_MA = canlab_results_fmridisplay(meta_fearcond_obj, 'compact', 'nooutline', 'trans')
drawnow, snapnow
clf()

%% display appetitive CSplusGTCSminus predictive weight maps 
%% sample 1
display_0188_svm_2 = orthviews(CSpGTCSm_0188_svm_obj);
drawnow, snapnow
clf()
display_obj_0188_svm = canlab_results_fmridisplay(CSpGTCSm_0188_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow

%% sample 2
display_0205_svm_2 = orthviews(CSpGTCSm_0205_svm_obj)
drawnow, snapnow
clf()
display_obj_0205_svm = canlab_results_fmridisplay(CSpGTCSm_0205_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow

%% sample 3
display_0187_svm_2 = orthviews(CSpGTCSm_0187_svm_obj)
drawnow, snapnow
clf()
display_obj_0187_svm = canlab_results_fmridisplay(CSpGTCSm_0187_svm_obj, 'compact', 'nooutline', 'trans');
drawnow, snapnow

%% display tables of region for all samples 
CSpGTCSm_0188_svm_obj_fdr_05 = threshold(CSpGTCSm_0188_svm_obj, 0.05, 'fdr');
CSpGTCSm_0205_svm_obj_fdr_05 = threshold(CSpGTCSm_0205_svm_obj, 0.05, 'fdr');
CSpGTCSm_0187_svm_obj_fdr_05 = threshold(CSpGTCSm_0187_svm_obj, 0.05, 'fdr');

%% Sample 1 regions
%% fdr .05 corr
r_CSpGTCSm_0188_svm_fdr_05 = region(CSpGTCSm_0188_svm_obj_fdr_05);    
    [rCSp0188_svm_fdr_05, rCSm0188_svm_fdr_05] = table(r_CSpGTCSm_0188_svm_fdr_05);       % add labels
drawnow, snapnow

%% unc .05
r_CSpGTCSm_0188_svm = region(CSpGTCSm_0188_svm_obj);    
    [rCSp0188_svm, rCSm0188_svm] = table(r_CSpGTCSm_0188_svm);       % add labels
drawnow, snapnow

%T = table(r_CSpGTCSm_0188_svm);
%writetable(T, fullfile(base_dir, '{scripts/html/r_CSpGTCSm_0188_svm.xlsx'))
%% Sample 2 regions
%% fdr .05 corr
r_CSpGTCSm_0205_svm_fdr_05 = region(CSpGTCSm_0205_svm_obj_fdr_05);    
    [rCSp0205_svm_fdr_05, rCSm0205_svm_fdr_05] = table(r_CSpGTCSm_0205_svm_fdr_05);       % add labels
drawnow, snapnow

%% unc .05
r_CSpGTCSm_0205_svm = region(CSpGTCSm_0205_svm_obj);    
    [rCSp0205_svm, rCSm0205_svm] = table(r_CSpGTCSm_0205_svm);       % add labels
drawnow, snapnow

%% Sample 3 regions
%% fdr .05 corr
r_CSpGTCSm_0187_svm_fdr_05 = region(CSpGTCSm_0187_svm_obj_fdr_05);    
    [rCSp0187_svm_fdr_05, rCSm0187_svm_fdr_05] = table(r_CSpGTCSm_0187_svm_fdr_05);       % add labels
drawnow, snapnow

%% unc .05
r_CSpGTCSm_0187_svm = region(CSpGTCSm_0187_svm_obj);    
    [rCSp0187_svm, rCSm0187_svm] = table(r_CSpGTCSm_0187_svm);       % add labels
drawnow, snapnow

%% contrast plots
% ------------------------------------------------------------------------
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
    disp('colormap range of blobs in case legend gets fucked up ')
    disp(o2.activation_maps{1}.cmaprange)
    drawnow, snapnow
    
end

