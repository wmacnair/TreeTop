%% treetop_layout: do layouts for this run
function [] = treetop_layout(input_struct, options_struct)
	% check inputs
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);
	
	% do force-directed graph layout
	calc_force_directed_layout(input_struct, options_struct);

	% which markers show greatest differences between branches?
	calc_most_significant_branching_markers(input_struct, options_struct)
end

%% calc_force_directed_layout: 
function calc_force_directed_layout(input_struct, options_struct)
	fprintf('calculating TreeTop layouts\n')

	% unpack
	output_dir 		= input_struct.output_dir;
	save_stem 		= input_struct.save_stem;

	% get graph structure
	graph_struct	= get_graph_struct(input_struct);
	n_graphs 		= size(graph_struct, 2);

	% define storage variable
	G_cell 			= cell(n_graphs, 1);
	layout_cell 	= cell(n_graphs, 1);
	title_cell 		= cell(n_graphs, 1);

	if options_struct.pool_flag
		% set up reproducible rng
		spmd
			cmrg 		= RandStream('mrg32k3a', 'seed', options_struct.seed);
			RandStream.setGlobalStream(cmrg);
		end

		% loop through all graphs
		parfor ii = 1:n_graphs
			% set up random stream
			s 				= RandStream.getGlobalStream();
			s.Substream 	= ii;

			% define graph to use,calculate layout
			G				= graph_struct(ii).inv_adj_matrix;
			layout 			= harel_koren_layout_faster(G);

			% store outputs
			G_cell{ii} 		= G;
			layout_cell{ii} = layout;
			title_cell{ii} 	= graph_struct(ii).name;
		end
	else
		% loop through all graphs
		for ii = 1:n_graphs
			% set up random stream
			rng(ii);

			% define graph to use,calculate layout
			G				= graph_struct(ii).inv_adj_matrix;
			layout 			= harel_koren_layout_faster(G);

			% store outputs
			G_cell{ii} 		= G;
			layout_cell{ii} = layout;
			title_cell{ii} 	= graph_struct(ii).name;
		end
	end

	% assemble into one figure
	put_all_figures_together(G_cell, layout_cell, title_cell, input_struct)
end

%% do_layout: 
function do_layout(G, layout, plot_title)
	% prepare plot outputs
	[ii jj ss] 			= find(G);
	[mm nn] 			= size(G);
	adjacency_matrix 	= sparse(ii, jj, repmat(1, numel(ii), 1), mm, nn);

	% plot
	hold on
	grey_gplot(G, layout);
	plot(layout(:, 1), layout(:, 2), '.', 'markersize', 20);
	hold off
	xlabel('TreeTop 1')
	ylabel('TreeTop 2')
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	title(plot_title, 'interpreter', 'none');
end

%% grey_gplot: 
function h2 = grey_gplot(layout_graph, layout_xy)
	gplot(layout_graph, layout_xy, '-k');
	h 				= gca;
	h2 				= get(h, 'Children');
	grey_val 		= 0.8;
	set(h2, 'color', [grey_val, grey_val, grey_val]);
end

%% put_all_figures_together: 
function put_all_figures_together(G_cell, layout_cell, title_cell, input_struct)
	% unpack
	output_dir 		= input_struct.output_dir;
	save_stem 		= input_struct.save_stem;

	% set up figure
	fig 			= figure('visible', 'off');
	n_rows 			= 2;
	n_cols 			= 3;

	% plot each layout
	for ii = 1:numel(G_cell)
		% do this one
		subplot(n_rows, n_cols, ii);
		do_layout(G_cell{ii}, layout_cell{ii}, title_cell{ii})
	end

	% save result as png
	name_stem 		= fullfile(output_dir, sprintf('%s all layouts', save_stem));
	fig_size 		= [12 8];
	plot_fig(fig, name_stem, 'png', fig_size)
end

%% calc_most_significant_branching_markers: 
function calc_most_significant_branching_markers(input_struct, options_struct)
	% get treetop outputs
	treetop_struct 		= get_treetop_outputs(input_struct);

	% calculate mean branching distance from this node to all other nodes
	mean_tree_dist 		= calculate_dists_from_branching_point(input_struct, options_struct, treetop_struct);

	% do ANOVA on the induced branches
	anova_results 		= identify_de_markers_with_anova(input_struct, options_struct, treetop_struct, mean_tree_dist);

	% save outputs
	save_branching_markers(mean_tree_dist, anova_results, input_struct)
end

%% calculate_dists_from_branching_point: calculates mean mst distance matrix for branching point
function [mean_tree_dist] = calculate_dists_from_branching_point(input_struct, options_struct, treetop_struct)
	% get tree distance bits
	tree_cell 	= treetop_struct.tree_cell;
	best_branches 	= treetop_struct.best_branches;
	branch_point 	= find(best_branches == 0);
	if length(branch_point) ~= 1
		error('something wrong with branches')
	end
	
	% do dijkstra for each tree
	n_nodes 		= options_struct.n_ref_cells;
	n_trees 		= options_struct.n_trees;

	fprintf('calculating mean mst distances from branching point\n');
	% define output array
	all_dijk_dists 	= zeros(n_trees, n_nodes);

	% loop
	for ii 	= 1:n_trees
		all_dijk_dists(ii, :) 	= dijkstra(tree_cell{ii}, branch_point);
	end

	% take mean
	mean_tree_dist 	= mean(all_dijk_dists, 1);

	% scale to have max distance 1
	max_dist 		= max(mean_tree_dist);
	mean_tree_dist 	= mean_tree_dist / max_dist;

	% have biggest branch on LHS (i.e. negative distances)
	one_idx 					= best_branches == 1;
	mean_tree_dist(one_idx) 	= -mean_tree_dist(one_idx);
end

%% identify_de_markers_with_anova: find which markers are most strongly differentially expressed between branches
function [anova_results] = identify_de_markers_with_anova(input_struct, options_struct, treetop_struct, mean_tree_dist)
	fprintf('calculating markers which are significantly different between branches\n')
	% unpack
	best_branches 		= treetop_struct.best_branches;
	used_values 		= treetop_struct.used_values;
	extra_values 		= treetop_struct.extra_values;
	n_used 				= size(used_values, 2);
	n_extra 			= size(extra_values, 2);

	% check we can actually do this
	group_size_check 	= check_anova_sizes(best_branches);
	if ~group_size_check
		% skip
		anova_results = struct( ...
			'anova_used', 	ones(1, n_used), ...
			'anova_extra', 	ones(1, n_extra) ...
			);
		fprintf('at least some branches too small to calculate which markers differentially expressed on branches\n')
		return
	end

	% do ANOVA to find which markers show greatest difference between branches
	branch_idx 		= best_branches ~= 0;
	anova_used 		= do_one_anova(branch_idx, best_branches, used_values);
	if ~isempty(extra_values)
		anova_extra 	= do_one_anova(branch_idx, best_branches, extra_values);
	else
		anova_extra 	 = [];
	end

	% store outputs
	anova_results = struct( ...
		'anova_used', 	anova_used, ...
		'anova_extra', 	anova_extra ...
		);
end

%% check_anova_sizes: 
function [group_size_check] = check_anova_sizes(best_branches)
	% check that it is ok to do anova (need n >= # groups + 1)
	[branch_count, labels] 	= grpstats(best_branches, best_branches, {'numel', 'gname'});
	labels 					= cellfun(@str2num, labels);

	% remove branch point
	keep_idx 				= labels ~= 0;
	branch_count 			= branch_count(keep_idx);
	labels 					= labels(keep_idx);

	% check number of branches against branch sizes
	n_branches 				= numel(labels);
	group_size_check 		= all( branch_count > n_branches + 1 );
end

%% do_one_anova: 
function [p_vals] = do_one_anova(branch_idx, best_branches, val_matrix)
	branch_vals 	= val_matrix(branch_idx, :);
	n_sel			= size(val_matrix, 2);
	actual_branches = best_branches(branch_idx);
	p_vals 			= arrayfun( @(ii) anova1(branch_vals(:, ii), actual_branches, 'off'), 1:n_sel );

	% correct for multiple testing
	p_vals 			= p_vals / size(val_matrix, 2);
end

%% save_branching_markers: 
function save_branching_markers(mean_tree_dist, anova_results, input_struct)
	output_file 	= fullfile(input_struct.output_dir, sprintf('%s marker significance.mat', input_struct.save_stem));
	save(output_file, 'mean_tree_dist', 'anova_results')
end
