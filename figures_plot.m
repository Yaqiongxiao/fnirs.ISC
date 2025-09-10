%% -------------------- plot interaction  --------------------
% Design: 2 (Cond Type: Recorded, Live) x 3 (Valence: Positive, Negative, Neutral)
% Both factors are within-subjects
clear;
clc;
close all;
datapath = 'xxx';
var_name = 'isc_persubject';
Cond_groups = {
    struct('name', 'Mon-social context', 'files', {{'ISC_recorded_pos_58sub.mat', 'ISC_recorded_neg_58sub.mat', 'ISC_recorded_neu_58sub.mat'}}), ...
    struct('name', 'Social context', 'files', {{'ISC_live_pos_58sub.mat', 'ISC_live_neg_58sub.mat', 'ISC_live_neu_58sub.mat'}})
    };
valence_cond_names = {'Positive', 'Negative', 'Neutral'};
Condtion_names = {'Non-social context', 'Social context'};

data_record_pos = load(fullfile(datapath, Cond_groups{1}.files{1})).(var_name)(:);
data_record_neg = load(fullfile(datapath, Cond_groups{1}.files{2})).(var_name)(:);
data_record_neu = load(fullfile(datapath, Cond_groups{1}.files{3})).(var_name)(:);

data_live_pos = load(fullfile(datapath, Cond_groups{2}.files{1})).(var_name)(:);
data_live_neg = load(fullfile(datapath, Cond_groups{2}.files{2})).(var_name)(:);
data_live_neu = load(fullfile(datapath, Cond_groups{2}.files{3})).(var_name)(:);

n_subj = length(data_record_pos);
if ~all(cellfun(@length, {data_record_neg, data_record_neu, data_live_pos, data_live_neg, data_live_neu}) == n_subj)
    error('The number of participants across all conditions must be equal！');
end
fprintf('Data loaded successfully. Total participants: %d.\n', n_subj);
SubjID = (1:n_subj)';
data_table = table(SubjID, ...
    data_record_pos, data_record_neg, data_record_neu, ...
    data_live_pos, data_live_neg, data_live_neu);

data_table.Properties.VariableNames = {'SubjID', ...
    'Record_Positive', 'Record_Negative', 'Record_Neutral', ...
    'Live_Positive', 'Live_Negative', 'Live_Neutral'};

figure('Color', 'w', 'Position', [200, 200, 800, 600], 'Name', '2x3 ANOVA');
means = [mean(data_table.Record_Positive), mean(data_table.Record_Negative), mean(data_table.Record_Neutral);
    mean(data_table.Live_Positive),  mean(data_table.Live_Negative),  mean(data_table.Live_Neutral)];

ses = [std(data_table.Record_Positive)/sqrt(n_subj), std(data_table.Record_Negative)/sqrt(n_subj), std(data_table.Record_Neutral)/sqrt(n_subj);
    std(data_table.Live_Positive)/sqrt(n_subj),  std(data_table.Live_Negative)/sqrt(n_subj),  std(data_table.Live_Neutral)/sqrt(n_subj)];
hold on;
errorbar(1:3, means(1,:), ses(1,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'CapSize', 10, 'DisplayName', Condtion_names{1});
errorbar(1:3, means(2,:), ses(2,:), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'CapSize', 10, 'DisplayName', Condtion_names{2});
hold off;
box off;
grid off;
xlim([0.5, 3.5]);
xticks(1:3);
xticklabels(valence_cond_names);
ylabel('ISC','FontName','Arial','FontWeight', 'bold','FontSize', 18);
legend('boxoff','show', 'Location', 'best', 'FontSize', 14);
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 16;
ax.LineWidth = 1.5;
title('Cond × Valence','FontName','Arial','FontSize', 20, 'FontWeight', 'bold');


%% -------------------- plot permutations  --------------------
clear;
clc;
close all;
input_folder = 'xxx';
save_path = 'xxx';
plot_save_path = fullfile(save_path, 'Permutation_Plots');
if ~exist(plot_save_path, 'dir')
    mkdir(plot_save_path);
end
groups = {'recorded', 'live'};
subject_counts = [58, 58];
conditions = {'neg', 'neu', 'pos'};
n_permutations = 1000;
valence_cond_names = {'Negative', 'Neutral', 'Positive'};
Condtion_names = {'Non-social', 'Social'};
for i = 1:length(groups)
    current_group = groups{i};
    nsub = subject_counts(i);
    if strcmp(current_group, 'recorded')
        display_group_name = Condtion_names{1}; % Non-social context
    else
        display_group_name = Condtion_names{2}; % 'live'，Social context
    end
    
    for j = 1:length(conditions)
        current_condition = conditions{j};
        display_valence_name = valence_cond_names{j};
        fprintf('\n--- Processing: Group [%s], Condition [%s] ---\n', current_group, current_condition);
        data_filename = sprintf('NULL_ISC_%s_%s_%dsub.mat', current_group, current_condition, nsub);
        full_data_path = fullfile(save_path, data_filename);
        
        fprintf('  Attempting to load: %s\n', full_data_path);
        if ~exist(full_data_path, 'file')
            warning('File not found, skipping: %s', full_data_path);
            continue;
        end
        
        loaded_data = load(full_data_path, 'isc_true_per_subject', 'null_isc_distribution');
        isc_true_per_subject = loaded_data.isc_true_per_subject;
        null_isc_distribution = loaded_data.null_isc_distribution;
        mean_isc_true = mean(isc_true_per_subject);
        mean_null_dist = mean(null_isc_distribution, 1);
        
        p_value = sum(mean_null_dist >= mean_isc_true) / n_permutations;
        
        h_fig = figure('Color', 'w', 'Name', sprintf('%s - %s Permutation Test', display_group_name, display_valence_name));
        hold on;
        hist(mean_null_dist, 50);
        h_patch = findobj(gca, 'Type', 'patch');
        if ~isempty(h_patch)
            set(h_patch, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
        end
        current_ylim = ylim;
        line([mean_isc_true, mean_isc_true], current_ylim, 'Color', 'r', 'LineWidth', 2.5);
        hold off;
        box off;
        grid off;
        xlabel('Mean ISC', 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Frequency', 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
        
        if p_value < 0.001
            p_value_str = 'p < 0.001';
        else
            p_value_str = sprintf('p = %.3f', p_value);
        end
        
        title({sprintf('%s - %s', display_group_name, display_valence_name), ...
            p_value_str}, ...
            'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
        
        legend('off');
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 14;
        ax.FontName = 'Arial';
        
        %         % --- save figures ---
        %         plot_output_filename = sprintf('Permutation_Plot_%s_%s.png', current_group, current_condition);
        %         full_plot_save_file = fullfile(plot_save_path, plot_output_filename);
        %         fprintf('  Saving plot to: %s\n', full_plot_save_file);
        %         saveas(h_fig, full_plot_save_file);
        %         close(h_fig);
    end
end


%%  -------------------- plot partial correlation  --------------------
clear;
clc;
close all;
datapath = 'xxx';
behavioral_file = 'ISC_58subj_info.csv';
var_name = 'isc_persubject';
control_vars_names = {'age', 'sexN', 'HI'};
record_files = {'ISC_recorded_pos_58sub.mat', 'ISC_recorded_neg_58sub.mat', 'ISC_recorded_neu_58sub.mat'};
live_files = {'ISC_live_pos_58sub.mat', 'ISC_live_neg_58sub.mat', 'ISC_live_neu_58sub.mat'};

data_record_pos = load(fullfile(datapath, record_files{1})).(var_name)(:);
data_record_neg = load(fullfile(datapath, record_files{2})).(var_name)(:);
data_record_neu = load(fullfile(datapath, record_files{3})).(var_name)(:);
data_live_pos = load(fullfile(datapath, live_files{1})).(var_name)(:);
data_live_neg = load(fullfile(datapath, live_files{2})).(var_name)(:);
data_live_neu = load(fullfile(datapath, live_files{3})).(var_name)(:);

isc_record_avg = mean([data_record_pos, data_record_neg, data_record_neu], 2);
isc_live_avg = mean([data_live_pos, data_live_neg, data_live_neu], 2);
delta_isc = isc_record_avg - isc_live_avg;
T_behav = readtable(fullfile(datapath, behavioral_file), 'TreatAsMissing', 'NULL');
all_behav_vars_names = T_behav.Properties.VariableNames;
interest_vars_names = setdiff(all_behav_vars_names, [control_vars_names, {'ID'}], 'stable');
control_matrix = T_behav{:, control_vars_names};

for i = 1:length(interest_vars_names)
    
    current_behav_var_name = interest_vars_names{i};
    current_behav_var = T_behav.(current_behav_var_name);
    
    [r, p] = partialcorr(delta_isc, current_behav_var, control_matrix, 'Rows', 'complete');
    
    if p < 0.05
        fprintf('  ** Results are significant. Generating partial correlation plots **\n');
        data_for_plot = [delta_isc, current_behav_var, control_matrix];
        complete_cases = ~any(isnan(data_for_plot), 2);
        
        isc_clean = delta_isc(complete_cases);
        behav_clean = current_behav_var(complete_cases);
        controls_clean = control_matrix(complete_cases, :);
        
        lm_isc = fitlm(controls_clean, isc_clean);
        res_isc = lm_isc.Residuals.Raw;
        lm_behav = fitlm(controls_clean, behav_clean);
        res_behav = lm_behav.Residuals.Raw;
        
        figure('Color', 'w', 'Name', sprintf('Partial Corr: Story Delta vs %s', current_behav_var_name));
        ax = gca;
        hold on;
        scatter(res_behav, res_isc, 80, 'filled', 'MarkerFaceColor', [0.2, 0.5, 0.8], 'MarkerFaceAlpha', 0.7);
        mdl_res = fitlm(res_behav, res_isc);
        h_ci = plot(mdl_res);
        set(h_ci(1), 'Marker', 'none');
        set(h_ci(2), 'Color', [0.8, 0.2, 0.2], 'LineWidth', 2.5);
        set(h_ci(3), 'Color', [0.8, 0.2, 0.2], 'LineStyle', '--', 'LineWidth', 1.2);
        set(h_ci(4), 'Color', [0.8, 0.2, 0.2], 'LineStyle', '--', 'LineWidth', 1.2);
        legend('off');
        hold off;
        box off;
        ax.LineWidth = 1.5;
        ax.FontSize = 14;
        ax.FontName = 'Arial';
        
        annotation_str = sprintf('r = %.3f\np = %.3f', r, p);
        text(ax.XLim(1) + 0.05*range(ax.XLim), ax.YLim(2) - 0.1*range(ax.YLim), ...
            annotation_str, 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
        xlabel(sprintf('%s', strrep(current_behav_var_name, '_', ' ')), 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
        ylabel('Delta ISC', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
        title('');
    end
    
end