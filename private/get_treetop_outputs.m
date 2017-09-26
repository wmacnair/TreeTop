%% get_treetop_outputs: open up saved variables from trees section
function [treetop_struct] = get_treetop_outputs(input_struct)
	fprintf('getting treetop outputs\n');

	% unpack
	output_dir 			= input_struct.output_dir;
	save_stem 			= input_struct.save_stem;

	% get all data (= node_idx_array, tree_array, tree_dist_array, celltype_vector)
	tree_variables_file 	= fullfile(output_dir, 'tree_variables.mat');
	load(tree_variables_file, 'outlier_idx', 'downsample_idx', 'tree_cell', 'cell_assignments', 'tree_dist_cell', 'centroids_idx', 'node_idx_array', 'celltype_vector', 'density');

	% get weight values
	samples_file 		= fullfile(output_dir, sprintf('%s_sample_counts.txt', save_stem));
	samples_struct 		= importdata(samples_file);
	sample_counts 		= samples_struct.data;
	sample_names 		= samples_struct.colheaders;

	% get marker values
	used_file 			= fullfile(output_dir, sprintf('%s_mean_used_markers.txt', save_stem));
	used_struct 		= importdata(used_file);
	mean_used 			= used_struct.data;
	extra_file 			= fullfile(output_dir, sprintf('%s_mean_extra_markers.txt', save_stem));
	if exist(extra_file, 'file')
		extra_struct 		= importdata(extra_file);
		mean_extra 			= extra_struct.data;
	else
		mean_extra 			= [];
	end

	% get used marker values
	used_file 			= fullfile(output_dir, sprintf('%s_all_used_markers.mat', save_stem));
	load(used_file, 'used_data');
	% used_values 		= used_data(centroids_idx, :);
	
	% get extra marker values, if they exist
	extra_file 			= fullfile(output_dir, sprintf('%s_all_extra_markers.mat', save_stem));
	if exist(extra_file, 'file')
		load(extra_file, 'extra_data');
		% extra_data 		= extra_data(centroids_idx, :);
	else
		extra_data 		= [];
	end

	% convert celltype_vector to nominal
	celltype_vector		= categorical(celltype_vector);

	% assemble into struct
	treetop_struct = struct( ...
		'tree_cell', 		{tree_cell}, ...
		'tree_dist_cell', 	{tree_dist_cell}, ...
		'used_values', 		{mean_used}, ...
		'extra_values', 	{mean_extra}, ...
		'sample_counts', 	{sample_counts}, ...
		'sample_names', 	{sample_names}, ...
		'centroids_idx', 	{centroids_idx}, ...
		'outlier_idx', 		{outlier_idx}, ...
		'downsample_idx', 	{downsample_idx}, ...
		'cell_assignments', {cell_assignments}, ...
		'density', 			{density}, ...
		'node_idx_array', 	{node_idx_array}, ...
		'used_data', 		{used_data}, ...
		'extra_data', 		{extra_data}, ...
		'celltype_vector',	{celltype_vector} ...
		);

	% if done, get scores and branches
	scores_file 		= fullfile(output_dir, sprintf('%s branching scores.txt', save_stem));
	if exist(scores_file, 'file')
		scores_struct 		= importdata(scores_file);
		branch_scores 		= scores_struct.data;

		% get branches
		branches_file 		= fullfile(output_dir, sprintf('%s best branches.txt', save_stem));
		branches_struct 	= importdata(branches_file);
		best_branches 		= branches_struct.data;

		% add to treetop_struct
		treetop_struct.branch_scores 	= branch_scores;
		treetop_struct.best_branches 	= best_branches;
	end

	% if done, get significance calculations for each marker
	output_file 	= fullfile(input_struct.output_dir, sprintf('%s marker significance.mat', input_struct.save_stem));
	if exist(output_file, 'file')
		load(output_file, 'mean_tree_dist', 'anova_results')
		treetop_struct.mean_tree_dist 	= mean_tree_dist;
		treetop_struct.anova_used 		= anova_results.anova_used;
		treetop_struct.anova_extra 		= anova_results.anova_extra;
	end

	% check whether cell_assignments_top exists
	file_details 		= who('-file', tree_variables_file);
	if ismember('cell_assignments_top', file_details)
		load(tree_variables_file, 'cell_assignments_top')
		treetop_struct.cell_assignments_top 	= cell_assignments_top;
	end
	if ismember('branch_parent_point', file_details)
		load(tree_variables_file, 'branch_parent_point')
		treetop_struct.branch_parent_point 	= branch_parent_point;
	end
end
