%% treetop_recursive: First sample ensemble of trees, do layouts. Then check whether we think there are further bifurcations.
function treetop_recursive(input_struct, options_struct)
	fprintf('\nrunning TreeTop recursively\n')

	% check whether pool is present
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);
	pool_check(options_struct);

	% do recursion bit
	treetop_one_recursion(input_struct, options_struct)

	% draw results together, plot
	recursive_struct 	= process_recursive_outputs(input_struct, options_struct);

	% do some plotting
	plot_recursive_outputs(input_struct, options_struct, recursive_struct);

	fprintf('\ndone.\n\n')
end

%% treetop_one_recursion: 
function treetop_one_recursion(input_struct, options_struct)
	% check inputs
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);

	% define things that might already be save
	tree_vars_file 		= fullfile(input_struct.output_dir, 'tree_variables.mat');
	branching_file 		= fullfile(input_struct.output_dir, sprintf('%s branching scores.txt', input_struct.save_stem));
	if ~exist(tree_vars_file, 'file') | ~exist(branching_file, 'file')
		% sample ensemble of trees
		treetop_trees(input_struct, options_struct)

		% calculate bifurcation score
		treetop_branching_scores(input_struct, options_struct)
	end

	% get treetop outputs
	treetop_struct 	= get_treetop_outputs(input_struct);

	% check whether branch has score above non-branching distribution
	[has_branch, n_branches] 	= check_for_branches(input_struct, options_struct, treetop_struct);

	if has_branch
		fprintf('branch point found for %s; zooming in to %d branches\n', input_struct.save_stem, n_branches);
		% split and save
		[branch_input_cell, branch_options_cell, branch_check] 		= split_treetop_branches(input_struct, options_struct, treetop_struct);

		% tidy up some memory before we recurse
		clear treetop_struct

		for ii = 1:length(branch_input_cell)
			% run recursion on this one
			branch_input 	= branch_input_cell{ii};
			branch_options 	= branch_options_cell{ii};

			% run recursive on this
			if branch_check(ii)
				treetop_one_recursion(branch_input, branch_options);
			else
				save_dummy_outputs_for_tiny_branches(branch_input, branch_options);
			end
		end
	else
		fprintf('no branch point found for %s; stopping recursion\n', input_struct.save_stem)
	end
end

%% check_for_branches: 
function [has_branch, n_branches] = check_for_branches(input_struct, options_struct, treetop_struct)
	% get non-branching distribution
	[n_points, n_dims] 	= size(treetop_struct.used_data);
	n_ref_cells			= options_struct.n_ref_cells;
	non_branching_distn = get_non_branching_distn(n_ref_cells, n_points, n_dims);

	% get scores, normalize
	branch_scores 		= treetop_struct.branch_scores;
	cutoff_q 			= quantile(non_branching_distn, options_struct.p_cutoff);
	normed_scores 		= branch_scores / cutoff_q;
	best_score 			= max(normed_scores);
	
	% check whether higher than 1
	has_branch 			= best_score > 1;

	% calculate n_branches
	n_branches 			= numel(unique(treetop_struct.best_branches)) - 1 ;
end

%% split_treetop_branches: 
function [branch_input_cell, branch_options_cell, branch_check] = split_treetop_branches(input_struct, options_struct, treetop_struct)
	% get branches, scores
	best_branches 		= treetop_struct.best_branches;
	branch_point 		= find(best_branches == 0);
	unique_branches 	= setdiff(unique(best_branches), 0);
	n_branches 			= length(unique_branches);
	branch_input_cell 	= cell(n_branches, 1);
	branch_options_cell = cell(n_branches, 1);
	branch_check		= true(n_branches, 1);

	% unpack
	save_stem			= input_struct.save_stem;
	cell_assignments 	= treetop_struct.cell_assignments;
	used_data 			= treetop_struct.used_data;
	extra_data 			= treetop_struct.extra_data;
	celltype_vector 	= treetop_struct.celltype_vector;

	% we want to keep track of the labels at the top level for all cells
	if isfield(treetop_struct, 'cell_assignments_top')
		cell_assignments_top 	= treetop_struct.cell_assignments_top;
	else
		cell_assignments_top 	= cell_assignments;
	end

	% remove filenames etc from input_struct
	if isfield(input_struct, 'filenames')
		input_struct	= rmfield(input_struct, 'filenames');
		input_struct	= rmfield(input_struct, 'file_annot');
	end

	% cycle through branches
	for ii = 1:n_branches
		% restrict to this branch
		this_branch 			= unique_branches(ii);

		% which datapoints to keep? (we also pass the branch point through)
		node_idx 				= find(ismember(best_branches, [0, this_branch]));
		nodes_to_keep 			= ismember(cell_assignments, node_idx);
		
		% check whether there are enough
		n_nodes					= sum(nodes_to_keep);
		if n_nodes < 200
			branch_check(ii)		= false;
		end

		% make new input_struct
		branch_input 	= input_struct;
		if regexp(save_stem, '^branch_')
			branch_stem		= sprintf('%s_%d', save_stem, this_branch);
		else
			branch_stem		= sprintf('branch_%d', this_branch);
		end
		branch_input.save_stem 	= branch_stem;

		% make new directory
		branch_dir 				= fullfile(input_struct.output_dir, branch_stem);
		if ~exist(branch_dir, 'dir')
			mkdir(branch_dir)
		end
		branch_input.output_dir = branch_dir;
		
		% restrict to this branch
		branch_used 			= used_data(nodes_to_keep, :);
		if ~isempty(extra_data)
			branch_extra 			= extra_data(nodes_to_keep, :);
		else
			branch_extra 			= [];
		end
		branch_celltypes 		= celltype_vector(nodes_to_keep);
		branch_assign_top 		= cell_assignments_top(nodes_to_keep);

		% which ones correspond to the parent branch point?
		branch_parent_point 	= cell_assignments(nodes_to_keep) == branch_point;

		% put into struct
		all_struct = struct( ...
			'all_data', 			{[branch_used, branch_extra]}, ...
			'all_labels', 			{branch_celltypes}, ...
			'all_markers', 			{{input_struct.used_markers{:}, input_struct.extra_markers{:}}}, ...
			'cell_assignments_top', {branch_assign_top}, ...
			'branch_parent_point', 	{branch_parent_point} ...
			);

		% save required outputs in new directory
		branch_file 			= sprintf('%s treetop inputs.mat', branch_stem);
		branch_path 			= fullfile(branch_dir, branch_file);
		save(branch_path, 'all_struct');

		% add this location to branch_input
		branch_input.data_dir 	= branch_dir;
		branch_input.mat_file 	= branch_file;
		
		% define options for this branch, including whether to do this branch
		branch_options 			= make_branch_options(all_struct, options_struct);

		% store this branch_input
		branch_input_cell{ii} 	= branch_input;
		branch_options_cell{ii} = branch_options;
	end
end

%% make_branch_options: check # ref cells: 
function branch_options = make_branch_options(all_struct, options_struct)
	% unpack
	n_ref_cells 	= options_struct.n_ref_cells;
	n_total 		= size(all_struct.all_data, 1);

	% checks there are at least 10 cells per node
	ratio_threshold 	= 10;
	cell_ratio 			= n_total / n_ref_cells;
	ratio_check 		= cell_ratio >= ratio_threshold;
	if ~ratio_check
		% change to smallest multiple of 10 achieving this threshold
		n_ref_cells 	= floor(n_total / ratio_threshold / 10 ) * 10;
	end

	% check there are at least 20 nodes; use this to decide whether to do this branch or not
	n_ref_cutoff 	= 20;
	n_ref_check 	= n_ref_cells >= n_ref_cutoff;

	% make options
	branch_options 				= options_struct;
	branch_options.n_ref_cells 	= n_ref_cells;
	branch_options.n_ref_check 	= n_ref_check;
	branch_options.outlier 		= 0;
end

%% save_dummy_outputs_for_tiny_branches:
function save_dummy_outputs_for_tiny_branches(branch_input, branch_options)
	% define things that might already be saved
	tiny_branch_file 	= fullfile(branch_input.output_dir, 'outputs_for_tiny_branch.mat');

	% save if necessary
	% if ~exist(tiny_branch_file, 'file')
	if true
		% unpack
		branch_dir 		= branch_input.output_dir;
		branch_stem 	= branch_input.save_stem;

		% load up inputs
		inputs_file 			= sprintf('%s treetop inputs.mat', branch_stem);
		inputs_path 			= fullfile(branch_dir, inputs_file);
		load(inputs_path);

		% fake tree_variables file
		cell_assignments_top 	= all_struct.cell_assignments_top;
		cell_assignments 		= cell_assignments_top;

		% save
		save(tiny_branch_file, 'cell_assignments_top', 'cell_assignments');
	end
end

%% process_recursive_outputs: rejoin all recursive outputs, somehow...
function recursive_struct = process_recursive_outputs(input_struct, options_struct)
	% recursively look for branches, get outputs for each
	fprintf('getting outputs for each recursive branch\n')
	% if already done, load
	recursive_output_file 	= fullfile(input_struct.output_dir, sprintf('%s recursive outputs.mat', input_struct.save_stem));
	% if exist(recursive_output_file, 'file')
	if false
		load(recursive_output_file, 'recursive_struct')
	else
		[stem_outputs, struct_outputs] 	= get_all_branch_outputs(input_struct.output_dir);
		% edit top
		stem_outputs{1}		= 'top';

		% put everything together, do some double-checking along the way
		recursive_struct 					= assemble_struct_outputs(stem_outputs, struct_outputs, input_struct, options_struct);
		recursive_struct.celltype_vector 	= struct_outputs{1}.celltype_vector;
		recursive_struct.cell_assignments 	= struct_outputs{1}.cell_assignments;
	end
end

%% get_all_branch_outputs: 
function [stem_outputs, struct_outputs] = get_all_branch_outputs(input_dir)
	% what is there in this directory?
	dir_details 		= dir(input_dir);
	dir_details 		= dir_details([dir_details.isdir]);
	dir_names 			= {dir_details.name};
	
	% get the outputs from this level
	[~, this_stem, ext] 	= fileparts(input_dir);
	if ~isempty(ext)
		this_stem 	= [this_stem, ext];
	end
	input_struct 		= struct('output_dir', input_dir, 'save_stem', this_stem);
	this_struct 		= get_treetop_outputs_recursive(input_struct);

	% put outputs together
	stem_outputs 		= {this_stem};
	struct_outputs 		= {this_struct};

	% find those which start with 'branch'
	branch_boolean 		= ~cellfun(@isempty, regexp(dir_names, '^branch'));
	if sum(branch_boolean) > 1
		% if we have branch directories, then go further down
		branch_dirs 					= cellfun(@(str) fullfile(input_dir, str), dir_names(branch_boolean), 'unif', false);
		[next_outputs, next_structs] 	= cellfun(@(this_dir) get_all_branch_outputs(this_dir), branch_dirs, 'unif', false);
		% mess about with cells
		if numel(next_structs) > 1
			% next_structs				= cellfun( @(this_cell) this_cell{:}, next_structs, 'unif', false);
			next_structs				= [ next_structs{:} ];
		end

		% turn into contiguous cell outputs
		stem_outputs 	= [stem_outputs{:}, next_outputs{:}];
		struct_outputs 	= {struct_outputs{:}, next_structs{:}};
	end
end

%% get_treetop_outputs_recursive: 
function this_struct = get_treetop_outputs_recursive(input_struct)
	% define file to check
	tiny_branch_file 	= fullfile(input_struct.output_dir, 'outputs_for_tiny_branch.mat');

	% if the branch was too small to run, load up the dummy outputs
	if exist(tiny_branch_file, 'file')
		load(tiny_branch_file, 'cell_assignments_top', 'cell_assignments');
		this_struct 	= struct( ...
			'cell_assignments_top', {cell_assignments_top}, ...
			'cell_assignments', 	{cell_assignments} ...
			);
		this_struct.run_flag 	= false;

	% otherwise load real outputs
	else
		this_struct 			= get_treetop_outputs(input_struct);
		this_struct.run_flag 	= true;

	end
end

%% assemble_struct_outputs: 
function recursive_struct = assemble_struct_outputs(stem_outputs, struct_outputs, input_struct, options_struct)
	% unpack
	n_nodes				= options_struct.n_ref_cells;

	% calculate branch depth for each branch
	branch_depths 		= cellfun(@(str) length(regexp(str, '_')), stem_outputs);
	max_depth 			= max(branch_depths);

	% trim names
	short_branch_names 	= cellfun(@(str) regexprep(str, 'branch_', ''), stem_outputs, 'unif', false);

	if max_depth == 0
		% we don't need to do all the stuff below if there's no branching
		[final_labels, branch_tree]		= make_no_branch_variables(n_nodes);
	else
		% define storage
		all_branch_labels 	= cell(n_nodes, max_depth);

		% start our tree
		branch_tree 		= initialize_branch_tree(stem_outputs);

		% cycle through branch depths
		for ii = 1:max_depth
			% restrict to this depth
			depth_idx 			= branch_depths == ii;
			depth_structs 		= struct_outputs( depth_idx );
			depth_names 		= short_branch_names( depth_idx );

			% join all together
			top_assigns_cell 	= cellfun(@(this_struct) this_struct.cell_assignments_top, depth_structs', 'unif', false);
			top_assigns_names 	= arrayfun(@(ii) repmat(depth_names(ii), 1, size(top_assigns_cell{ii}, 1)), 1:length(top_assigns_cell), 'unif', false);
			top_assigns 		= cell2mat(top_assigns_cell);
			top_assigns_names 	= [top_assigns_names{:}]';

			% do crosstab to count
			[depth_branch_votes, ~, ~, labels] 	= crosstab(top_assigns, top_assigns_names);
			depth_labels 						= labels(:, 2);
			depth_labels 						= depth_labels(~cellfun(@isempty, depth_labels));

			% define labels
			[~, max_branch_idx] 				= max(depth_branch_votes, [], 2);
			node_labels 						= arrayfun(@(idx) depth_labels{idx}, max_branch_idx, 'unif', false);

			% store these in the right place
			node_idx 							= cellfun(@str2num, labels(:, 1));
			all_branch_labels(node_idx, ii) 	= node_labels;

			% % do some tweaking of branch_tree
			% branch_tree			= update_branch_tree(ii, branch_depths, short_branch_names, struct_outputs, branch_tree);
		end

		% make tree symmetric
		branch_tree 		= branch_tree + branch_tree';

		% fill in missing labels where necessary
		top_best_branches					= struct_outputs{1}.best_branches;
		missing_idx							= find(cellfun(@isempty, all_branch_labels(:, 1)));
		all_branch_labels(missing_idx, 1)	= arrayfun(@(idx) num2str(top_best_branches(idx)), missing_idx, 'unif', false);

		% label each cell according to deepest level
		empty_matrix 		= cellfun(@isempty, all_branch_labels);
		node_max_depth 		= arrayfun(@(ii) max(find(~empty_matrix(ii, :))), 1:n_nodes);
		final_labels 		= arrayfun(@(ii) all_branch_labels{ii, node_max_depth(ii)}, 1:n_nodes, 'unif', false);
	end

	% make branching point lookup table
	[branch_points, point_scores] 	= cellfun(@(this_struct) find_branch_point_xys(this_struct, options_struct), struct_outputs);
	
	% one option: 
		% for every cell
		% get deepest label possible
		% then label each node according to most popular

		% somehow also some checking?
			% have we covered all nodes?
			% are any inappropriately duplicated?
	% assemble output
	recursive_struct 	= struct( ...
		'node_labels',		{final_labels}, ...
		'branch_tree', 		{branch_tree}, ...
		'branch_points', 	{branch_points}, ...
		'branch_names', 	{short_branch_names}, ...
		'point_scores',		{point_scores} ...
		);

	% save outputs
	recursive_output_file 	= fullfile(input_struct.output_dir, sprintf('%s recursive outputs.mat', input_struct.save_stem));
	save(recursive_output_file, 'recursive_struct')
end

%% make_no_branch_variables: 
function [final_labels, branch_tree] = make_no_branch_variables(n_nodes)
	final_labels	= repmat({'1'}, 1, n_nodes);
	branch_tree		= 0;
end

%% initialize_branch_tree: outputs vector of parents for each node
function branch_tree = initialize_branch_tree(stem_outputs)
	% tweak first entry
	stem_outputs{1} 	= 'branch';

	% define storage
	n_branches 			= numel(stem_outputs);
	branch_tree 		= zeros(n_branches);

	for ii = 2:n_branches
		% get this branch
		this_branch			= stem_outputs{ii};

		% trim, match
		end_idx 			= regexp(this_branch, '_[0-9]+$')-1;
		trimmed_branch		= this_branch(1:end_idx);
		parent_idx			= find(strcmp(trimmed_branch, stem_outputs));

		% store
		branch_tree(parent_idx, ii)	= 1;
	end
end

%% update_branch_tree: 
% if at top level, do nothing. otherwise, for each sub-branch:
% find which branch has biggest connection to previous branch point
% need branch flag for parent run (i.e. label for each cell,
% stating whether it's part of the parent branch point. then let
% the branches vote. this branch is connected to self and to 
function branch_tree = update_branch_tree(ii, branch_depths, short_branch_names, struct_outputs, branch_tree)
	% first allocate all branches to parent
	if ii > 1
		% find parent branch points for each one at this depth
		depth_idx 			= branch_depths == ii;
		% depth_names 		= short_branch_names(depth_idx);
		% parent_branches		= unique(cellfun(@(this_name) this_name(1:end-2), depth_names, 'unif', false));
		depth_names 		= short_branch_names(depth_idx);
		end_idxs 			= cellfun(@(this_name) regexp(this_name, '_[0-9]+$')-1, depth_names);
		parent_branches		= unique(cellfun(@(ii) depth_names{ii}(1:end_idxs(ii)), 1:length(depth_names), 'unif', false));

		for this_branch = parent_branches
			% which is parent, and which are children
			this_branch			= this_branch{1};
			parent_idx			= find(strcmp(this_branch, short_branch_names));

			% get parent outputs to decide split
			parent_struct 		= struct_outputs{parent_idx};
			
			% unpack a bit
			branch_parent_point = parent_struct.branch_parent_point;
			cell_assignments 	= parent_struct.cell_assignments;
			best_branches 		= parent_struct.best_branches;

			% exclude outliers
			outlier_idx 		= cell_assignments == 0;
			if sum(outlier_idx) > 0
				fprintf('we have outliers')
				branch_parent_point 	= branch_parent_point(~outlier_idx);
				cell_assignments 		= cell_assignments(~outlier_idx);
			end

			% count how much of the parent branch point ends up in each child branch
			branches_by_cell 							= best_branches(cell_assignments);
			[branch_point_votes, branch_name_checks] 	= grpstats(branch_parent_point, branches_by_cell, {'sum', 'gname'});
			if strcmp(branch_name_checks{1}, '0')
				branch_point_votes 	= branch_point_votes(2:end);
				branch_name_checks 	= branch_name_checks(2:end);
			else
				error('first branch name should be 0')
			end

			% find max
			[~, max_branch] 	= max(branch_point_votes);
			max_branch_name 	= sprintf('%s_%d', this_branch, max_branch);

			% find where this is in tree
			max_child_idx 		= find(strcmp(max_branch_name, short_branch_names));
			grandparent_idx 	= find(branch_tree(:, parent_idx));

			% do editing
			branch_tree(grandparent_idx, parent_idx) 	= 0;
			branch_tree(parent_idx, max_child_idx) 		= 0;
			branch_tree(grandparent_idx, max_child_idx)	= 1;
			branch_tree(max_child_idx, parent_idx)		= 1;
		end
	end
end

%% find_branch_point_xys: 
function [branch_point, point_score] = find_branch_point_xys(this_struct, options_struct)
	% check if treetop was run on this branch
	if this_struct.run_flag == true
		% unpack
		branch_scores 				= this_struct.branch_scores;
		cell_assignments 			= this_struct.cell_assignments;

		% get appropriate non-branching distribution, normalize scores
		[n_ref_cells, n_dims] 		= size(this_struct.used_values);
		n_points					= size(this_struct.cell_assignments, 1);
		non_branching_distn 		= get_non_branching_distn(n_ref_cells, n_points, n_dims);
		cutoff_q 					= quantile(non_branching_distn, options_struct.p_cutoff);
		normed_scores 				= branch_scores / cutoff_q;

		% now calculate location of best branch point
		if isfield(this_struct, 'cell_assignments_top')
			% unpack
			cell_assignments_top 		= this_struct.cell_assignments_top;

			% calculate where the branch point lines up at the top
			[max_score, max_idx] 		= max(normed_scores);
			branch_idx 					= cell_assignments == max_idx;
			branch_top 					= cell_assignments_top(branch_idx);

			% find which node is best fit
			[node_counts, node_labels] 	= grpstats(branch_top, branch_top, {'numel', 'gname'});
			[~, max_node] 				= max(node_counts);
			branch_point 				= str2num(node_labels{max_node});
		else
			[max_score, branch_point]	= max(normed_scores);
		end
		point_score 			= max_score;

    else
		% dummy branch_point
		branch_point 			= 0;

		% give branch score of 0
		point_score 			= 0;
	end
end

%% plot_recursive_outputs: 
function [] = plot_recursive_outputs(input_struct, options_struct, recursive_struct)
	fprintf('plotting outputs for recursive TreeTop\n')

	% get layout
	recursive_flag 	= true;
	layout_struct 	= get_layout_struct(input_struct, options_struct, recursive_flag);

	% unpack
	node_labels 	= recursive_struct.node_labels;
	branch_points 	= recursive_struct.branch_points;
	branch_names 	= recursive_struct.branch_names;
	point_scores 	= recursive_struct.point_scores;
	branch_tree 	= recursive_struct.branch_tree;

	celltype_vector 	= recursive_struct.celltype_vector;
	cell_assignments 	= recursive_struct.cell_assignments;
	n_samples 			= length(unique(celltype_vector));

	file_ext 		= options_struct.file_ext;

	layout_xy 		= layout_struct.layout_xy;
	layout_graph 	= layout_struct.layout_graph;

	% edit where branch points are for branches where treetop not run
	branch_points 	= tweak_branch_point_xys(layout_xy, node_labels, branch_names, branch_points);
	branch_xy 		= layout_xy(branch_points, :);

	% set up figure
	fig 			= figure('visible', 'off');
	n_cols 			= 2;
	n_rows 			= 2;
	plot_ii 		= 1;
	grey_val 		= 0.8;

	% plot branches as colour
	subplot(n_rows, n_cols, plot_ii)
	plot_ii 		= plot_ii+1;
	plot_branches_coloured(layout_graph, layout_xy, grey_val, node_labels)

	% plot branches as text
	subplot(n_rows, n_cols, plot_ii)
	plot_ii 		= plot_ii+1;
	plot_branches_labelled(layout_graph, layout_xy, grey_val, node_labels)

	% plot contingency table
	if n_samples > 1
		subplot(n_rows, n_cols, plot_ii)
		plot_ii 		= plot_ii+1;
		plot_contingency_table_recursive(recursive_struct)
	end

	% plot branch points
	subplot(n_rows, n_cols, plot_ii)
	plot_ii 		= plot_ii+1;
	plot_branch_tree(layout_graph, layout_xy, grey_val, branch_tree, branch_xy, point_scores)

	% save outputs
	plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s recursive branches', input_struct.save_stem));
	plot_unit 		= 4;
	fig_size 		= [plot_unit*n_cols*1.1, plot_unit*n_rows];
	plot_fig(fig, plot_stem, file_ext, fig_size)
end

%% tweak_branch_point_xys: 
function branch_points = tweak_branch_point_xys(layout_xy, node_labels, branch_names, branch_points)
	for ii = 1:length(branch_points)
		if branch_points(ii) == 0
			% find which points correspond to this branch
			this_name 			= branch_names{ii};
	
			% find point closest to centre
			branch_idx 			= find(strcmp(node_labels, this_name));
			branch_xy 			= layout_xy(branch_idx, :);
			mean_xy 			= mean(branch_xy, 1);
			temp 				= branch_xy - repmat(mean_xy, size(branch_xy, 1), 1);
			[~, near_idx] 		= min( sum(temp.^2,2) );
			branch_points(ii) 	= branch_idx(near_idx);
		end
	end
end

%% plot_branch_tree: 
function plot_branch_tree(layout_graph, layout_xy, grey_val, branch_tree, branch_xy, point_scores)
	hold on
	% plot whole graph behind
	gplot(layout_graph, layout_xy, '-k');
	h 			= gca;
	h2 			= get(h, 'Children');
	set(h2, 'color', [grey_val, grey_val, grey_val])

	% plot graph connecting branch points
	gplot(branch_tree, branch_xy, '-k');

	% plot branch scores in order
	[~, point_idx] 	= sort(point_scores);
	scatter(branch_xy(point_idx, 1), branch_xy(point_idx, 2), [], point_scores(point_idx)', 'filled');
	score_range 		= [0, ceil(max([point_scores(:); 1]))];
	% point_size		= 40;
	% size_vector 	= point_size * ones(size(branch_xy, 1), 1);
	% size_vector(1) 	= point_size * 2;
	plot_size 		= get(gca, 'Position');
	bar_obj 		= colorbar;
	set(gca, 'Position', plot_size);
	bar_pos 		= get(bar_obj, 'position');
	bar_pos(3:4) 	= bar_pos(3:4) / 2;
	bar_pos(1)	 	= bar_pos(1) - bar_pos(3);
	bar_pos(2)	 	= bar_pos(2) + bar_pos(4)/2;
	set(bar_obj, 'position', bar_pos)
	ylabel(bar_obj, 'Relative branching score', 'interpreter', 'none')

	% where the position arguments are [xposition yposition width height].
	caxis( score_range )
	xlim([0,1])
	ylim([0,1])
	xlabel('TreeTop 1'); ylabel('TreeTop 2')

	% labels
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	hold off
end

%% plot_branches_coloured: 
function plot_branches_coloured(layout_graph, layout_xy, grey_val, node_labels)
	hold on
	gplot(layout_graph, layout_xy, '-k');
	h 				= gca;
	h2 				= get(h, 'Children');
	grey_val 		= 0.8;
	set(h2, 'color', [grey_val, grey_val, grey_val])
	gscatter_for_recursive(layout_xy(:,1), layout_xy(:,2), node_labels);

	% % plot branch scores in order
	% [~, point_idx] 	= sort(point_scores);
	% gplot(branch_tree, branch_xy, '-k');
	% scatter(branch_xy(point_idx, 1), branch_xy(point_idx, 2), 60, point_scores(point_idx)', 'filled');
	% score_range 	= [0, ceil(max([point_scores(:); 1]))];

	% sort out plot
	xlim([0, 1])
	ylim([0, 1])
	xlabel('TreeTop 1'); ylabel('TreeTop 2')
	hold off

	% label graph
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
end

%% gscatter_for_recursive: slightly fancy plotting for recursive plots
function [h_legend] = gscatter_for_recursive(x, y, node_labels, options)
	% sort out holding
	hold_status 	= ishold;
	if ~hold_status
		hold on
	end

	% get palette
	ctype 		= 'qual';
	palette 	= 'Set1';
	point_size 	= 10;
	legend_flag = false;
	location 	= 'EastOutside';

	if nargin > 3 
		if isfield(options, 'palette')
			palette 	= options.palette;
		end
		if isfield(options, 'size')
			point_size 	= options.size;
		end
		if isfield(options, 'legend')
			legend_flag = options.legend;
		end
		if isfield(options, 'leg_pos')
			location 	= options.leg_pos;
		end
	end

	% regex to get top level of branches, then see where they are
	g_top		   = cellfun(@(c) regexp(c, '^[0-9]+', 'match'), node_labels);
	[top_vals, ~, g_idx] 	= unique(g_top);
	if ~iscell(top_vals)
		top_vals 		= arrayfun(@num2str, top_vals, 'unif', false);
	end

	% how many top branches, sub-branches?
	n_top 			= length(top_vals);
	n_all 			= length(unique(node_labels));

	% set up plots
	h 				= zeros(n_all, 1);

	% if more than nine of those, don't do anything fancy
	if n_top > 9
		[~, ~, all_idx] = unique(node_labels);
		scatter(x, y, point_size, g_idx, 'filled');
		
	else
		% define top-level palette
		n_pal 			= max(n_top, 3);
		pal_top 		= cbrewer(ctype, palette, n_pal);

		% define counter for top branches
		col_inc 		= 1;

		% loop through top branches
		for ii = 1:n_top
			% which branch is this?
			ii_val 			= top_vals{ii};
			branch_idx 		= g_idx == ii;

			% how many sub-branches?
			sub_vals 		= unique(node_labels(branch_idx));
			n_sub 			= length(sub_vals);

			if n_sub==1
				this_idx 		= strcmp(node_labels, sub_vals{1});
				h(col_inc) 		= plot(x(this_idx), y(this_idx), '.', 'color', pal_top(ii, :), 'markersize', point_size*2);
				col_inc 		= col_inc + 1;

			else 
				% make palette for sub_branches based on top branch colour
				hsv_top 		= rgb2hsv(pal_top(ii, :));
				hsv_sub 		= repmat(hsv_top, n_sub, 1);
				s_range 		= [0.95, 0.05];
				s_vector 		= linspace(s_range(1), s_range(2), n_sub);
				v_vector 		= 1 - 0.6*(1-max(s_range))./(1 - s_vector);
				hsv_sub(:, 2) 	= s_vector;
				hsv_sub(:, 3) 	= v_vector;
				pal_sub 		= hsv2rgb(hsv_sub);

				% plot each sub-branch individually
				for jj = 1:n_sub
					this_sub 		= sub_vals{jj};
					this_idx 		= strcmp(node_labels, this_sub);
					h(col_inc) 		= plot(x(this_idx), y(this_idx), '.', 'color', pal_sub(jj, :), 'markersize', point_size*2);
					col_inc 		= col_inc + 1;
				end
			end
		end
	end

	% add branch points on top?

	if legend_flag
		h_legend 		= legend(h, g_vals{:}, 'Location', location);
	else
		h_legend 		= [];
	end

	if ~hold_status
		hold off
	end
end

%% plot_branches_labelled: 
function plot_branches_labelled(layout_graph, layout_xy, grey_val, node_labels)
	hold on
	gplot(layout_graph, layout_xy, '-k');
	h 			= gca;
	h2 			= get(h, 'Children');
	set(h2, 'color', [grey_val, grey_val, grey_val])
	text(layout_xy(:,1), layout_xy(:,2), node_labels, 'fontsize', 6, 'interpreter', 'none', 'horizontalalignment', 'center', 'verticalalignment', 'middle');
	xlim([0, 1])
	ylim([0, 1])
	xlabel('TreeTop 1'); ylabel('TreeTop 2')
	hold off

	% label graph
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
end

%% plot_contingency_table_recursive: 
function plot_contingency_table_recursive(recursive_struct)
	% unpack
	cell_assignments 	= recursive_struct.cell_assignments;
	celltype_vector 	= recursive_struct.celltype_vector;
	node_labels 		= recursive_struct.node_labels;
	[branch_names, ~, node_idx] = unique(node_labels);

	% remove any outliers
	non_outlier_idx 	= cell_assignments ~= 0;
	cell_assignments 	= cell_assignments(non_outlier_idx);
	celltype_vector 	= celltype_vector(non_outlier_idx);

	% % assign branches to all original cells
	% branches_by_cell 	= node_idx(cell_assignments);

	% % recalculate counts
	% [branch_counts, count_labels] 	= grpstats(branches_by_cell, branches_by_cell, {'numel', 'gname'});
	% count_labels					= cellfun(@(str) str2num(str), count_labels);

	% % lookup table which has new label for each branch, according to order
	% [~, count_order] 				= sort(-branch_counts);
	% [~, ranked_labels] 				= sort(count_order);
	% branches_by_size				= branch_names(count_order);

	% % put in order from top left to top right, plot table
	% [mean_branch_by_celltype, labels]		= grpstats(ranked_labels(branches_by_cell), celltype_vector, {'mean', 'gname'});
	% [~, sort_idx]							= sort(mean_branch_by_celltype);
	% col_order								= labels(sort_idx);
	% plot_contingency_table(branch_names(branches_by_cell), celltype_vector, branches_by_size, col_order)

	% recalculate counts
	[branch_counts, count_labels] 	= grpstats(node_idx, node_idx, {'numel', 'gname'});
	count_labels					= cellfun(@(str) str2num(str), count_labels);

	% lookup table which has new label for each branch, according to order
	[~, count_order] 				= sort(-branch_counts);
	[~, ranked_labels] 				= sort(count_order);
	branches_by_size				= branch_names(count_order);

	% put in order from top left to top right, plot table
	[mean_branch_by_celltype, labels]		= grpstats(ranked_labels(node_idx(cell_assignments)), celltype_vector, {'mean', 'gname'});
	[~, sort_idx]							= sort(mean_branch_by_celltype);
	col_order								= labels(sort_idx);
	plot_contingency_table(branch_names(node_idx(cell_assignments)), celltype_vector, branches_by_size, col_order)
end