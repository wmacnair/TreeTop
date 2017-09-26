%% treetop_trees: (1) load data, calculate density, take subsample (2) identify outlier / downsample points
% (3) do kmeans++ seeding, allocate each point to closest seed (4) repeatedly: sample one point from each cluster, fit tree
% between them, record marker values and which celltypes these were (5) do layouts
function treetop_trees(input_struct, options_struct)
	% do seed
	rng(options_struct.seed);

	% load up data, calculate density if necessary, take subsample to make it more manageable
	[sample_struct, options_struct] = load_and_process_data(input_struct, options_struct);

	% find appropriate outlier and threshold points
	[outlier_idx, downsample_idx] 	= calc_outlier_downsample_idx(sample_struct, options_struct);

	% do k-means ++ seeding (start with small # of clusters, for ease of HK layouts)
	% also allocate each point to one of these clusters (basically like )
	[centroids_idx, cell_assignments] = partition_cells(sample_struct, outlier_idx, downsample_idx, options_struct);
	% error('stop here')

	% set up tree variables
	n_nodes 		= numel(centroids_idx);
	n_trees 		= options_struct.n_trees;
	tree_cell 		= cell(n_trees, 1);
	tree_dist_cell 	= cell(n_trees, 1);
	node_idx_array 	= NaN(n_trees, n_nodes);

	% set up reproducible rng
	if options_struct.pool_flag
		spmd
			cmrg 		= RandStream('mrg32k3a', 'seed', options_struct.seed);
			RandStream.setGlobalStream(cmrg);
		end

		% sample all the trees
		fprintf('sampling %d trees (. = 50):  ', n_trees)
		parfor ii = 1:n_trees
			% set up random stream
			s 			= RandStream.getGlobalStream();
			s.Substream = ii;
			if mod(ii, 50) == 0
				fprintf([char(8), '. '])
			end

			% sample one point from each cluster, connect as MST
			[this_node_idx, this_tree, this_tree_dist] 	= get_a_tree(sample_struct, cell_assignments, options_struct);
			
			% store output
			node_idx_array(ii, :) 	= this_node_idx;
			tree_cell{ii} 			= this_tree;
			tree_dist_cell{ii} 		= this_tree_dist;
		end
	else
		% sample all the trees
		fprintf('sampling %d trees (. = 50): ', n_trees)
		for ii = 1:n_trees
			% set up random stream
			rng(ii);

			if mod(ii, 50) == 0
				fprintf('.')
			end

			% sample one point from each cluster, connect as MST
			[this_node_idx, this_tree, this_tree_dist] 	= get_a_tree(sample_struct, cell_assignments, options_struct);
			
			% store output
			node_idx_array(ii, :) 	= this_node_idx;
			tree_cell{ii} 			= this_tree;
			tree_dist_cell{ii} 		= this_tree_dist;
		end
		fprintf('\n')
	end
	% save output as union / freq graph
	[union_graph, freq_union_graph] 	= calculate_union_graphs(tree_cell, tree_dist_cell, input_struct, options_struct);

	% save assignment of every cell to a reference cell
	save_cell_assignments(input_struct, sample_struct, cell_assignments);

	% save marker values
	save_marker_values(input_struct, sample_struct, centroids_idx, cell_assignments);

	% save tree values
	save_sampled_trees(input_struct, options_struct, outlier_idx, downsample_idx, centroids_idx, cell_assignments, node_idx_array, tree_cell, tree_dist_cell, sample_struct);
end

%% load_and_process_data: 
function [sample_struct, options_struct] = load_and_process_data(input_struct, options_struct)
	% open all files, stitch together (with sample labels)
	all_struct 		= get_all_files(input_struct);

	% possibly remove ball around zero
	all_struct 		= remove_zero_ball(all_struct, options_struct);

	% take subsample to make things more manageable (maybe 100k?)
	sample_struct 	= take_sample(all_struct, options_struct);

	% check # ref cells
	options_struct 	= check_ref_cells(all_struct, options_struct);

	% get density values
	sample_struct 	= add_density_values(sample_struct, options_struct);
end

%% take_sample: 
function sample_struct = take_sample(all_struct, options_struct)
	% skip if not specified
	if ~isfield(options_struct, 'sample_size')
		sample_struct 	= all_struct;
		return
	end

	% otherwise, unpack
	fprintf('taking smaller sample\n')
	sample_size 	= options_struct.sample_size;
	used_data 		= all_struct.used_data;
	extra_data 		= all_struct.extra_data;

	% pick sample
	n_cells			= size(used_data, 1);
	if sample_size < n_cells
		sample_idx 		= randsample(n_cells, sample_size);
	else
		fprintf('requested sample size (%d) >= # cells (%d); all cells taken as sample\n', sample_size, n_cells);
		sample_idx 		= 1:n_cells;
	end

	% restrict to this
	sample_used 				= used_data(sample_idx, :);
	sample_labels 				= all_struct.all_labels(sample_idx);
	if isempty(extra_data)
		sample_extra			= [];
	else
		sample_extra 			= extra_data(sample_idx, :);
	end

	% assemble output
	sample_struct 				= all_struct;
	sample_struct.used_data 	= sample_used;
	sample_struct.extra_data 	= sample_extra;
	sample_struct.all_labels 	= sample_labels;

	% do bit for keeping track of recursion
	if isfield(all_struct, 'cell_assignments_top')
		sample_struct.cell_assignments_top 		= all_struct.cell_assignments_top(sample_idx);
	end
	% do bit for keeping track of recursion
	if isfield(all_struct, 'branch_parent_point')
		sample_struct.branch_parent_point 		= all_struct.branch_parent_point(sample_idx);
	end
end

%% check_ref_cells: check # ref cells
function options_struct = check_ref_cells(all_struct, options_struct)
	% unpack
	n_ref_cells 	= options_struct.n_ref_cells;
	n_total 		= size(all_struct.used_data, 1);

	% checks there are at least 10 cells per node
	ratio_threshold = 10;
	cell_ratio 		= n_total / n_ref_cells;
	ratio_check 	= cell_ratio >= ratio_threshold;
	if ~ratio_check
		% change to smallest multiple of 10 achieving this threshold
		n_ref_cells 	= floor(n_total / ratio_threshold / 10 )*10;
	end

	% check there are at least 20 nodes
	n_ref_cutoff 	= 20;
	n_ref_check 	= n_ref_cells >= n_ref_cutoff;

	if ~n_ref_check
		if ratio_check
			fprintf('attempting to learn presence of branch points with only 20 nodes; stopping\n');
			fprintf('try increasing number of reference nodes\n');
		else
			fprintf('this seems to be too small a dataset for TreeTop\n');
		end
		error('too few reference cells')
	else
		if ~ratio_check
			fprintf('too few cells per ref node; reducing # ref nodes to %d\n', n_ref_cells);
		end
	end
end

%% add_density_values: add density values to the structure
function sample_struct = add_density_values(sample_struct, options_struct)
	% if no downsampling specified, don't do it
	if options_struct.outlier == 0 & options_struct.threshold == 1
		density_vector 		= ones(size(sample_struct.used_data, 1), 1);
	else
		% otherwise calculate density
		density_vector 			= calc_density(sample_struct, options_struct);
	end

	% add to sample_struct
	sample_struct.density 	= density_vector;
end

%% calc_outlier_downsample_idx:
function [outlier_idx, downsample_idx] = calc_outlier_downsample_idx(sample_struct, options_struct);
	if isempty(options_struct.outlier)
		options_struct.outlier 		= 0;
	end
	if isempty(options_struct.threshold)
		options_struct.threshold 	= 1;
	end

	% unpack
	density_vector 	= sample_struct.density;

	% exclude outlier
	outlier_idx 	= exclude_outliers(density_vector, options_struct);

	% downsample by density
	downsample_idx 	= downsample_to_target(density_vector, options_struct, outlier_idx);
end

%% exclude_outliers: exclude outlier
function [outlier_idx] 		= exclude_outliers(density_vector, options_struct)
	% unpack
	outlier 	= options_struct.outlier;

	% define any outliers to remove
	if outlier == 0
		% if 0, no outliers
		outlier_idx 	= [];
	else
		% find appropriate density value
		outlier_q 		= quantile(density_vector, outlier);
		outlier_idx 	= find(outlier_q > density_vector);
	end
end

%% downsample_to_target: downsample by density
function [downsample_idx] 	= downsample_to_target(density_vector, options_struct, outlier_idx)
	% unpack
	threshold 	= options_struct.threshold;

	% define any cells to downsample on basis of high density
	if threshold == 1
		% if 1, don't downsample
		downsample_idx 	= outlier_idx;

	else
		% histogram(density_vector)
		% find quantile
		threshold_q 	= quantile(density_vector, threshold);

		% define probabilities for keeping during downsampling
		keep_prob 		= min(1, threshold_q./density_vector);

		% exclude outliers
		keep_prob(outlier_idx) 	= 0;
		% so outlier_idx is always subset of downsample_idx

		% pick some cells
		rand_vector 	= rand(numel(density_vector), 1);
		downsample_idx 	= find(rand_vector > keep_prob);
	end
end

%% get_a_tree: 
function [this_node_idx, this_tree, this_tree_dist] = get_a_tree(sample_struct, cell_assignments, options_struct)
	% sample one point from each cluster
	cluster_list 	= setdiff(unique(cell_assignments), 0);
	this_node_idx 	= arrayfun(@(ii) sample_fn(cell_assignments, ii), cluster_list);

	% now get the corresponding locations in space
	cluster_data 	= sample_struct.used_data(this_node_idx, :);

	% and do a tree on it
	switch options_struct.metric_name
		case 'L1'
			[this_tree, this_tree_dist] = mst_expanded(cluster_data, 'L1');
		case 'L2'
			[this_tree, this_tree_dist] = mst_expanded(cluster_data, 'euclidean');
		case 'correlation'
			[this_tree, this_tree_dist] = mst_expanded(cluster_data, 'corr');
		case 'angle'
			[this_tree, this_tree_dist] = mst_expanded(cluster_data, 'angle');
		otherwise
			error('options_struct.metric_name is not a valid option')
	end
end

%% sample_fn: 
function sample_idx = sample_fn(cell_assignments, ii)
	cluster_idx 	= find(cell_assignments == ii);
	n_in_cluster 	= length(cluster_idx);

	sample_idx 		= cluster_idx(randsample(n_in_cluster, 1));
end

%% calculate_union_tree: superimpose all trees to calculate edge frequencies and mean edge lengths
function [union_graph, freq_union_graph] = calculate_union_graphs(tree_cell, tree_dist_cell, input_struct, options_struct)
	fprintf('calculating summaries of ensemble of trees\n')

	% unpack
	n_trees 				= options_struct.n_trees;
	n_ref_cells 			= options_struct.n_ref_cells;

	% define variables
	sparse_total_length 	= sparse(n_ref_cells, n_ref_cells);
	sparse_total_freq 		= sparse(n_ref_cells, n_ref_cells);

	% cycle through each cluster file
	for ii = 1:n_trees
		% get current one
		this_tree_dist 		= tree_dist_cell{ii};

		% update total length, total freq
		sparse_total_length = sparse_total_length + this_tree_dist;
		sparse_total_freq 	= sparse_total_freq + (this_tree_dist > 0);
	end

	% calculate sparse_mean_output, remove any NaNs resulting from division by zero
	union_graph				= sparse_total_length ./ sparse_total_freq;
	union_graph(find(isnan(union_graph))) = 0;
	
	% use frequency to define alternative distance
	[non_zero_ii, non_zero_jj] 	= find(sparse_total_freq);
	total_freq_vals 			= 1 ./ nonzeros(sparse_total_freq);
	freq_union_graph 			= sparse(non_zero_ii, non_zero_jj, total_freq_vals);

	% save both graphs
	tree_filename 				= fullfile(input_struct.output_dir, sprintf('%s_union_tree.mat', input_struct.save_stem));
	save(tree_filename, 'union_graph');
	tree_filename 				= fullfile(input_struct.output_dir, sprintf('%s_freq_union_tree.mat', input_struct.save_stem));
	save(tree_filename, 'freq_union_graph');
end

%% save_cell_assignments: save assignment of every cell to a reference cell
function [] = save_cell_assignments(input_struct, sample_struct, cell_assignments)
	fprintf('counting celltype composition of each reference cell\n')

	% convert to indexing variable
	[this_tab, ~, ~, this_labels] 	= crosstab(cell_assignments, sample_struct.all_labels);

	% outliers labelled as 0; remove these
	if strcmp(this_labels(1), '0')
		this_tab 	= this_tab(2:end, :);
	end
	non_empty						= find(cellfun(@(ii) ~isempty(ii), this_labels(:, 2)));
	this_header 					= this_labels(non_empty, 2)';

	% save output
	samples_file 	= fullfile(input_struct.output_dir, sprintf('%s_sample_counts.txt', input_struct.save_stem));
	save_txt_file(samples_file, this_header, this_tab);
end

%% save_marker_values: 
function [] = save_marker_values(input_struct, sample_struct, centroids_idx, cell_assignments)
	% unpack
	output_dir 		= input_struct.output_dir;
	save_stem 		= input_struct.save_stem;

	% restrict to cluster centres
	used_data 		= sample_struct.used_data;
	used_markers 	= sample_struct.used_markers;
	cluster_used 	= sample_struct.used_data(centroids_idx, :);
	extra_data 		= sample_struct.extra_data;

	% save outputs for used markers
	markers_file 	= fullfile(output_dir, sprintf('%s_mean_used_markers.txt', save_stem));
	save_txt_file(markers_file, sample_struct.used_markers, cluster_used);
	markers_file 	= fullfile(output_dir, sprintf('%s_all_used_markers.mat', save_stem));
	save(markers_file, 'used_markers', 'used_data');

	if ~isempty(extra_data)
		extra_markers 	= sample_struct.extra_markers;
		cluster_extra 	= sample_struct.extra_data(centroids_idx, :);
		% save extra marker values
		markers_file 	= fullfile(output_dir, sprintf('%s_mean_extra_markers.txt', save_stem));
		save_txt_file(markers_file, sample_struct.extra_markers, cluster_extra);
		markers_file 	= fullfile(output_dir, sprintf('%s_all_extra_markers.mat', save_stem));
		save(markers_file, 'extra_markers', 'extra_data');
	end
end

%% save_sampled_trees: 
function [] = save_sampled_trees(input_struct, options_struct, outlier_idx, downsample_idx, centroids_idx, cell_assignments, node_idx_array, tree_cell, tree_dist_cell, sample_struct)
	% restrict to just samples for label
	celltype_vector 	= sample_struct.all_labels;
	density 			= sample_struct.density;

	% save union graph as mat file
	tree_filename 		= fullfile(input_struct.output_dir, 'tree_variables.mat');

	% define which variables to save
	save_vars 			= {'outlier_idx', 'downsample_idx', 'density', 'centroids_idx', 'cell_assignments', 'node_idx_array', 'tree_cell', 'tree_dist_cell', 'celltype_vector', 'input_struct', 'options_struct'};
	if isfield(sample_struct, 'cell_assignments_top')
		cell_assignments_top	= sample_struct.cell_assignments_top;
		save_vars 				= {save_vars{:}, 'cell_assignments_top'};
	end
	if isfield(sample_struct, 'branch_parent_point')
		branch_parent_point		= sample_struct.branch_parent_point;
		save_vars 				= {save_vars{:}, 'branch_parent_point'};
	end

	% do saving
	save(tree_filename, save_vars{:})
end

