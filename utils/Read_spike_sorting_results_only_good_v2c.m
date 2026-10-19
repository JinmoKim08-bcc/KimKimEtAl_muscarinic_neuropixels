function [Data] = Read_spike_sorting_results_only_good_v2c(fname, Data, file_idx, param)

kilo_result = struct();
kilo_result.spike_clusters   = readNPY(fullfile(fname, 'spike_clusters.npy'));    % [nSpikes×1]
kilo_result.spike_amp        = readNPY(fullfile(fname, 'amplitudes.npy'));
kilo_result.spike_locs       = readNPY(fullfile(fname, 'spike_times.npy'));
kilo_result.cluster_info     = tdfread(fullfile(fname, 'cluster_info.tsv'));
kilo_result.cluster_group     = tdfread(fullfile(fname, 'cluster_group.tsv'));
kilo_result.templates        = readNPY(fullfile(fname, 'templates.npy'));        % [nTemplates×T×C]
kilo_result.spike_templates  = readNPY(fullfile(fname, 'spike_templates.npy'));  % [nSpikes×1]

fs = 30000;  % Hz

all_clusters = unique(kilo_result.spike_clusters);

%%% eliminate clusters labeled as noise
noise_group = kilo_result.cluster_group.group;
noise_clus_list = zeros(size(noise_group, 1), 1);
for clus_idx = 1:size(noise_group,1)
    if contains(string(noise_group(clus_idx, :)), 'noise')
        noise_clus_list(clus_idx) = 0; % it is noise
    elseif contains(string(noise_group(clus_idx, :)), 'mua')
        noise_clus_list(clus_idx) = 0; % it is multi unit
    else
        noise_clus_list(clus_idx) = 1; % it is (postulated) single unit
    end
end

%%% eliminate clusters which encounters ISI violations
thresh = 0.5;
isi_viol = kilo_result.cluster_info.isi_viol;
for clus_idx = 1:size(noise_group,1)
    if isi_viol(clus_idx) >= thresh
        noise_clus_list(clus_idx) = 0;
    end
end

%%% eliminate clusters which encounters amplitude_cutoff
amplitude_cutoff = 0.1;
amp_cut = kilo_result.cluster_info.amplitude_cutoff;
for clus_idx = 1:size(noise_group,1)
    if amp_cut(clus_idx) > amplitude_cutoff
        noise_clus_list(clus_idx) = 0;
    end
end

noise_clus_list = logical(noise_clus_list);
clus_id_list = 1:size(kilo_result.cluster_group.cluster_id, 1);
selected_cluster_id = clus_id_list(noise_clus_list);
kilo_result.selected_cluster_id = selected_cluster_id;
Data(file_idx).kilo_result = kilo_result;

% kilo_result.cluster_info.cluster_id(selected_cluster_id(141))

%%% Extract SU activity
SU = struct();
num_clus = numel(kilo_result.selected_cluster_id);
num_orig_clus = numel(clus_id_list);
for clus_idx = 1:num_clus
    SU(clus_idx).ch = kilo_result.cluster_info.ch(selected_cluster_id(clus_idx));
    SU(clus_idx).clus_id = kilo_result.cluster_group.cluster_id(selected_cluster_id(clus_idx));
    SU(clus_idx).locs = double(Data(file_idx).kilo_result.spike_locs(Data(file_idx).kilo_result.spike_clusters == SU(clus_idx).clus_id));
    SU(clus_idx).time = double(Data(file_idx).kilo_result.spike_locs(Data(file_idx).kilo_result.spike_clusters == SU(clus_idx).clus_id))/param.fs;
    SU(clus_idx).kilo_result = kilo_result;
end

Data(file_idx).num_clus = num_clus;
Data(file_idx).num_orig_clus = num_orig_clus;
Data(file_idx).SU = SU;

end

