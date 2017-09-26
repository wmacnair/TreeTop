%% partition_cells: given a list of points to be downsampled and a list of ouliers, identifies 
% centroids via k-means ++ seeding amongst those that have been downsampled (downsample_idx includes outlier_idx)
function [centroids_idx, cell_assignments] = partition_cells(sample_struct, outlier_idx, downsample_idx, options_struct)
	fprintf('selecting reference nodes\n')

	% unpack
	used_data 			= sample_struct.used_data;
	n_cells 			= size(used_data, 1);

	% make list of nodes to keep, select these
	keep_idx 			= setdiff(1:n_cells, downsample_idx);
	selected_data 		= used_data(keep_idx, :);
	n_selected 			= size(selected_data, 1);

	% find some start points
	kk 					= options_struct.n_ref_cells;
	seed_options 		= struct('metric', options_struct.metric_name, 'pool_flag', options_struct.pool_flag, 'verbose', false);
	[~, cents_idx_idx] 	= kmeans_plus_plus(selected_data, kk, kk, 3, seed_options, options_struct);
	cents_idx_idx 		= sort(cents_idx_idx);

	% turn this into indices for whole set rather than just downsampled cells
	centroids_idx 		= keep_idx(cents_idx_idx);

	% which points are not centroids?
	not_outliers 		= setdiff(1:n_cells, outlier_idx);
	not_centroids 		= setdiff(not_outliers, centroids_idx);

	% check that this makes sense
	partition_check 	= isequal(sort([outlier_idx(:); not_centroids(:); centroids_idx(:)])', 1:n_cells);
	if ~partition_check
		error('something went wrong in calculation of reference nodes')
	end

	% assign clusters to each
	X 					= used_data(centroids_idx, :);
	Y 					= used_data(not_centroids, :);

	% get distance matrix: # centroids * # datapoints
	fprintf('calculating distance from all cells to reference nodes\n')
	if options_struct.pool_flag
		D 				= all_distance_fn_par(X, Y, options_struct.metric_name);
	else
		D 				= all_distance_fn(X, Y, options_struct.metric_name);
	end
	fprintf('labelling each cell with closest reference node\n')
	cell_assignments 	= calc_cell_assignments(options_struct, D, centroids_idx, not_centroids, outlier_idx);
end

%% plot_centroid_marker_distns: 
function [] = plot_centroid_marker_distns(centroids_idx, sample_struct, paramset)
	% define parameters
	n_col 				= 4;
	edge_vector 		= -2:1:9;

	% unpack
	used_data 			= sample_struct.used_data;
	used_markers 		= sample_struct.used_markers;
	n_markers 			= size(used_data, 2);
	n_row 				= ceil(n_markers / n_col);

	% restrict to just these values
	centroid_data 		= used_data(centroids_idx, :);

	% plot histogram of each
	figure('name', 'Centroid univariate distributions')
	for ii = 1:n_markers
		subplot(n_row, n_col, ii);
		histogram(centroid_data(:, ii), edge_vector);
		title(used_markers{ii});
	end

	% 
	figure('name', 'Centroid bivariate distributions')
	% ii is the y axis marker
	for ii = 1:(n_markers-1)
		% jj is the y axis marker
		for jj = (ii+1):n_markers
			% which plot?
			subplot_idx = (ii-1)*(n_markers-1) + (jj-1);
			subplot(n_markers-1, n_markers-1, subplot_idx);
			% plot distn
			plot(centroid_data(:, jj), centroid_data(:, ii), '.');
			xlim([min(edge_vector), max(edge_vector)]);
			ylim([min(edge_vector), max(edge_vector)]);
			
			% label appropriately
			title_str 	= [used_markers{ii} ' vs ' used_markers{jj}];
			title(title_str);
		set(gca,'FontSize',8)
		end
	end

	% save outputs
	plot_file 	= fullfile(paramset.output_dir, 'marker biaxials.png');
	set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 20 20]);
	r = 300; % pixels per inch
	print(gcf, '-dpng', sprintf('-r%d', r), plot_file);
end

%% get_closest_cluster: finds which centroid each point is closest to
function [closest_cluster] = get_closest_cluster(this_column)
	% close_options 	= find( this_column == min(this_column) );
	% if numel(close_options) > 1
	% 	fprintf('joint closest!\n');
	% end
	% n_options 		= length(close_options);
	% selected_option = randsample(n_options, 1);
	% closest_cluster = close_options(selected_option);
	[~, closest_cluster] 	= min(this_column);
end

%% calc_cell_assignments: 
function [cell_assignments] = calc_cell_assignments(options_struct, D, centroids_idx, not_centroids, outlier_idx)
	% which is closest?
	n_non_centroids 	= numel(not_centroids);
	closest_clusters 	= NaN(n_non_centroids, 1);

	% loop!
	if options_struct.pool_flag
		parfor ii = 1:n_non_centroids
			closest_clusters(ii) 	= get_closest_cluster(D(:, ii));
		end
	else
		for ii = 1:n_non_centroids
			closest_clusters(ii) 	= get_closest_cluster(D(:, ii));
		end
	end

	% how many cells overall?
	n_cells 			= sum([numel(centroids_idx), numel(not_centroids), numel(outlier_idx)]);

	% put together into cell_assignments
	cell_assignments 					= NaN(n_cells, 1);
	cell_assignments(centroids_idx) 	= 1:length(centroids_idx);
	cell_assignments(not_centroids) 	= closest_clusters;
	cell_assignments(outlier_idx) 	= 0;
end
