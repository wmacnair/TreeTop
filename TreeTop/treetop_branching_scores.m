%% treetop_branching_scores: calculate bifurcation score for given treetop run
function [] = treetop_branching_scores(input_struct, options_struct)
	% check inputs
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);
	
	% define set of thresholds
	threshold_list 			= 0.99:-0.01:0.01;

	% get outputs we need
	treetop_struct 			= get_treetop_outputs(input_struct);

	% get consensus matrices for the branch points
	consistency_array 		= get_consistency_matrices(input_struct, options_struct, treetop_struct);
	
	% what happens when we gradually increase the threshold at which we keep edges?
	consensus_branch_sizes 	= calc_consensus_branch_sizes(consistency_array, input_struct, options_struct, threshold_list);

	% calculate mean size of these
	branching_scores 		= calc_branching_scores(consensus_branch_sizes, input_struct);

	% save best branches
	calc_and_save_best_branches(branching_scores, consistency_array, consensus_branch_sizes, threshold_list, input_struct, options_struct)
end

%% get_consistency_matrices: cut each tree at each node, take average across all trees
function consistency_array = get_consistency_matrices(input_struct, options_struct, treetop_struct, selected_nodes)

	% unpack
	n_nodes				= options_struct.n_ref_cells;
	n_trees				= options_struct.n_trees;
	tree_dist_cell		= treetop_struct.tree_dist_cell;

	% initialize various things (third dimension in arrays is cut point)
	consistency_array 	= zeros(n_nodes, n_nodes, n_nodes);

	% do bifurcation scores for each tree
	fprintf('splitting each tree at each node (%d trees, . = 50):  ', n_trees);
	if options_struct.pool_flag
		parfor ii = 1:n_trees
			if mod(ii, 50) == 0
				fprintf([char(8), '. '])
			end

			% get tree, calculate consensus matrix
			this_dist_tree 		= tree_dist_cell{ii};
			this_consistency	= calc_one_consistency_matrix(this_dist_tree);

			% store
			consistency_array 	= consistency_array + this_consistency;
		end
	else
		for ii = 1:n_trees
			if mod(ii, 50) == 0, 	fprintf('.'); 	end

			% get tree, calculate consensus matrix
			this_dist_tree 		= tree_dist_cell{ii};
			this_consistency	= calc_one_consistency_matrix(this_dist_tree);

			% store
			consistency_array 	= consistency_array + this_consistency;
		end
		fprintf('\n')
	end

	% normalize
	consistency_array 		= consistency_array / n_trees;

	% put -1s in columns and rows
	for ii = 1:n_nodes
		neg_idx 							= 1:n_nodes ~= ii;
		consistency_array(neg_idx, ii, ii) 	= -1;
		consistency_array(ii, neg_idx, ii) 	= -1;
	end
end

%% calc_one_consistency_matrix: 
function this_consistency = calc_one_consistency_matrix(this_dist_tree)
	% unpack, initialize
	n_nodes 			= size(this_dist_tree, 1);
	this_consistency		= zeros(n_nodes, n_nodes, n_nodes);
	this_counts			= zeros(n_nodes, n_nodes, n_nodes);

	% loop through all selected nodes
	for this_node = 1:n_nodes
		% make tree with this node removed
		split_tree 					= this_dist_tree;
		split_tree(this_node, :) 	= 0;
		split_tree(:, this_node) 	= 0;

		% calculate components
		blocks 						= components(split_tree);

		% update consensus matrix
		branch_list 				= unique(blocks);
		for kk = branch_list
			kk_idx 											= blocks == kk;
			this_consistency(kk_idx, kk_idx, this_node) 	= 1;
		end

		% but remove info from the node we cut at
		this_consistency(this_node, this_node, this_node) 	= 0;
	end
end

%% calc_consensus_branch_sizes: 
function [consensus_branch_sizes] = calc_consensus_branch_sizes(consistency_array, input_struct, options_struct, threshold_list)
	% unpack
	n_nodes 		= options_struct.n_ref_cells;

	% check if size file already exists
	size_file 		= fullfile(input_struct.output_dir, sprintf('%s consensus branch sizes.mat', input_struct.save_stem));

	% define thresholds, cutoffs to use
	threshold_max 	= 1;
	cutoff_list 	= threshold_max - threshold_list;
	n_thresholds 	= length(threshold_list);
	n_cutoffs 		= length(cutoff_list);

	% define storage variable
	consensus_branch_sizes 		= cell(n_nodes, 1);

	% do clustering etc for each consensus matrix
	fprintf('calculating branches with hierarchical clustering\n');
	if options_struct.pool_flag
		parfor this_node = 1:n_nodes
			% get relevant point and matrix
			this_consistency 	= squeeze(consistency_array(:, :, this_node));

			% remove negative values
			this_idx 			= 1:n_nodes == this_node;
			this_consistency 	= this_consistency( ~this_idx, ~this_idx );

			% convert similarity matrix into dissimilarity, then vector
			dis_mat 			= 1 - this_consistency;
			dis_vec 			= squareform(dis_mat, 'tovector');

			% apply linkage (i.e. group nodes according to linkage specified in link_str)
			this_link 			= linkage(dis_vec, 'single');

			% extract clusters at each cutpoint
			cluster_mat 		= cluster(this_link, 'cutoff', cutoff_list, 'criterion', 'distance');
			% how many clusters is this?
			this_max 			= max(cluster_mat(:));

			% count how many clusters at each cutoff, reorder by size
			[mesh_ii, mesh_jj] 	= meshgrid(1:n_cutoffs, 1:this_max);
			cluster_sizes 		= arrayfun( @(ii, jj) sum(cluster_mat(:, ii) == jj), mesh_ii, mesh_jj);
			cluster_sizes 		= cell2mat(arrayfun(@(ii) sort(cluster_sizes(:, ii), 'descend'), 1:n_cutoffs, 'unif', false))';
			
			% remove singletons, then remove any empty columns
			cluster_sizes(cluster_sizes == 1)	= 0;
			empty_cols							= sum(cluster_sizes, 1) == 0;
			cluster_sizes(:, empty_cols)		= [];

			% store
			consensus_branch_sizes{this_node} 	= cluster_sizes;
		end
	else
		for this_node = 1:n_nodes
			% get relevant point and matrix
			this_consistency 	= squeeze(consistency_array(:, :, this_node));

			% remove negative values
			this_idx 			= 1:n_nodes == this_node;
			this_consistency 	= this_consistency( ~this_idx, ~this_idx );

			% convert similarity matrix into dissimilarity, then vector
			dis_mat 			= 1 - this_consistency;
			dis_vec 			= squareform(dis_mat, 'tovector');

			% apply linkage (i.e. group nodes according to linkage specified in link_str)
			this_link 			= linkage(dis_vec, 'single');

			% extract clusters at each cutpoint
			cluster_mat 		= cluster(this_link, 'cutoff', cutoff_list, 'criterion', 'distance');
			% how many clusters is this?
			this_max 			= max(cluster_mat(:));

			% count how many clusters at each cutoff, reorder by size
			[mesh_ii, mesh_jj] 	= meshgrid(1:n_cutoffs, 1:this_max);
			cluster_sizes 		= arrayfun( @(ii, jj) sum(cluster_mat(:, ii) == jj), mesh_ii, mesh_jj);
			cluster_sizes 		= cell2mat(arrayfun(@(ii) sort(cluster_sizes(:, ii), 'descend'), 1:n_cutoffs, 'unif', false))';
			
			% remove singletons, then remove any empty columns
			cluster_sizes(cluster_sizes == 1)	= 0;
			empty_cols							= sum(cluster_sizes, 1) == 0;
			cluster_sizes(:, empty_cols)		= [];

			% store
			consensus_branch_sizes{this_node} 	= cluster_sizes;
		end
	end

	% save outputs
	save(size_file, 'consensus_branch_sizes', 'threshold_list', 'threshold_max');
end

%% calc_branching_scores: 
function branching_scores = calc_branching_scores(consensus_branch_sizes, input_struct)
	% unpack, initialize
	n_nodes 			= numel(consensus_branch_sizes);
	branching_scores 	= zeros(n_nodes, 1);

	% exclude clusters accounting for 1% or less of nodes
	cutoff 				= 0.01;

	% loop through all
	fprintf('calculating branching scores\n');
	for ii = 1:n_nodes
		% get this branch size value
		this_branch_sizes 	= consensus_branch_sizes{ii};
		this_branch_sizes 	= this_branch_sizes / n_nodes * 100;

		% check whether there are at least 3 clusters
		n_clusters 			= size(this_branch_sizes, 2);

		% if not, set to 0
		if n_clusters <= 2
			branching_scores(ii) 	= 0;

		% use mean size of all clusters above third, above a given threshold
		else
			% apply cutoff
			this_branch_sizes( this_branch_sizes(:)<=cutoff ) 	= 0;

			% define cols to use
			smaller_branches		= sum(this_branch_sizes(:, 3:end), 2);
			branching_scores(ii)	= mean(smaller_branches);
		end
	end

	% save bifurcation scores
	fprintf('saving branching scores\n');
	scores_file 	= fullfile(input_struct.output_dir, sprintf('%s branching scores.txt', input_struct.save_stem));
	save_txt_file(scores_file, {'score'}, branching_scores);
end

%% calc_and_save_best_branches: finds best branches, saves them
function calc_and_save_best_branches(branching_scores, consistency_array, consensus_branch_sizes, threshold_list, input_struct, options_struct)
	% unpack
	n_nodes 			= options_struct.n_ref_cells;

	% get best one
	[~, best_point]		= max(branching_scores);
	best_point 			= best_point(1);
	best_boolean 		= 1:n_nodes == best_point;

	% pick threshold giving largest 3rd branch
	this_sizes 			= consensus_branch_sizes{best_point};
	if size(this_sizes, 2) < 3
		best_branches		= ones(size(branching_scores));
		fprintf('no branches found at all for run %s\n', input_struct.save_stem);
	else
		[~, size_idx] 		= max( this_sizes(:,3) );
		best_thresh 		= threshold_list(size_idx);

		% get and process consensus matrix
		this_consistency 	= squeeze(consistency_array(:, :, best_boolean));
		this_consistency 	= this_consistency( ~best_boolean, ~best_boolean );

		% do clustering
		best_branches 		= cluster_consistency_matrix(this_consistency, best_thresh);

		% put in sensible order which also excludes singletons
		best_branches 		= process_best_branches(best_branches, this_consistency, best_boolean, options_struct);
	end

	% save
	branches_file 		= fullfile(input_struct.output_dir, sprintf('%s best branches.txt', input_struct.save_stem));
	save_txt_file(branches_file, {'branch'}, best_branches)
end

%% cluster_consistency_matrix: 
function best_branches = cluster_consistency_matrix(this_consistency, best_thresh)
	% convert similarity matrix into dissimilarity, then vector
	dis_mat 			= 1 - this_consistency;
	dis_vec 			= squareform(dis_mat, 'tovector');

	% apply linkage (i.e. group nodes via single linkage)
	this_link 			= linkage(dis_vec, 'single');
	best_branches 		= cluster(this_link, 'cutoff', 1 - best_thresh, 'criterion', 'distance');
end

%% process_best_branches: we want to remove singletons from our clustering
function [best_branches] = process_best_branches(best_branches, this_consistency, best_boolean, options_struct)
	% get components, and their sizes
	[counts, labels] 	= grpstats(best_branches, best_branches, {'numel', 'gname'});
	labels 				= cellfun(@(str) str2num(str), labels);

	% option to force 3 branches to be returned
	if options_struct.three_flag
		% find top 3 clusters
		[~, count_idx] 	= sort(-counts);
		n_counts 		= length(counts);
		top_3 			= min(n_counts, 3);
		multiples		= labels(count_idx(1:top_3));
	else
		% find non-singleton labels
		cutoff			= 0.01;
		multiples_idx	= counts / sum(counts) > cutoff & counts > 1;
		multiples		= labels(multiples_idx);
	end
	
	% get node locations
	singles				= setdiff(labels, multiples);
	singles_idx			= find(ismember(best_branches, singles));
	multiples_idx		= find(ismember(best_branches, multiples));

	% find nearest multiple for every single
	simil_matrix			= this_consistency( singles_idx, multiples_idx );
	[~, closest_multiple]	= max(simil_matrix, [], 2);

	% relabel
	branches_no_singles					= best_branches;
	branches_no_singles(singles_idx)	= best_branches(multiples_idx(closest_multiple));
	% relabel to shorter list
	[~, ~, branches_no_singles]			= unique(branches_no_singles);

	% recalculate counts
	[new_counts, new_labels] 	= grpstats(branches_no_singles, branches_no_singles, {'numel', 'gname'});
	new_labels					= cellfun(@(str) str2num(str), new_labels);

	% lookup table which has new label for each branch, according to order
	[~, count_order] 			= sort(-new_counts);
	[~, ranked_labels] 			= sort(count_order);
	branches_by_size			= ranked_labels(branches_no_singles);

	% check it has worked
	check_table		= tabulate(branches_by_size);
	if any( check_table(2:end, 1) - check_table(1:end-1, 1) ~= 1)
		error('relabelling went wrong');
	end
	if any( check_table(2:end, 2) > check_table(1:end-1, 2) )
		error('relabelling went wrong');
	end
	
	% put together into full branch, with 0 at branch point
	best_branches 					= zeros(size(best_boolean))';
	best_branches(best_boolean) 	= 0;
	best_branches(~best_boolean) 	= branches_by_size;
end
