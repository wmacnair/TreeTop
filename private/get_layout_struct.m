%% get_layout_struct: 
function layout_struct = get_layout_struct(input_struct, options_struct, layout_tree_idx, recursive_flag)
	% define default input
	if ~exist('recursive_flag', 'var')
		recursive_flag 	= false;
	end

	% double check inputs
	[input_struct, options_struct] 		= check_treetop_inputs(input_struct, options_struct);

	% define save file
	layout_file 		= fullfile(input_struct.output_dir, sprintf('best layout %s.mat', input_struct.save_stem));

	% if idx given, pick which type of layout to use, check input is ok, and force computation of layout
	if exist('layout_tree_idx', 'var') && ~isempty(layout_tree_idx)
		if ~ismember(layout_tree_idx, 1:6)
			error('invalid layout_tree_idx')
		end
	% otherwise use previously calculated layout if it exists
	else
		layout_tree_idx 	= 6;

		% check if exists
		if exist(layout_file, 'file')
			% if exists, check the date
			tree_vars_file 	= fullfile(input_struct.output_dir, 'tree_variables.mat');
			tree_vars_info 	= dir(tree_vars_file);
			layout_info 	= dir(layout_file);

			if tree_vars_info.datenum < layout_info.datenum | recursive_flag
				load(layout_file, 'layout_struct');
				return
			end
		end
	end


	% get graphs
	graph_struct 		= get_graph_struct(input_struct);

	% get details
	layout_graph_name 	= graph_struct(layout_tree_idx).name;
	layout_graph 		= graph_struct(layout_tree_idx).inv_adj_matrix;

	% run thing
	fprintf('doing layout for %s\n', layout_graph_name);
	if isfield(options_struct, 'layout_seed')
		layout_seed 		= options_struct.layout_seed;
	else
		layout_seed 		= 1;
	end
	layout_xy 			= get_best_layout(layout_graph, layout_seed, options_struct);

	layout_struct 	= struct( ...
		'layout_xy', 	{layout_xy}, ...
		'layout_graph', {layout_graph} ...
		);

	save(layout_file, 'layout_struct');
end

%% get_best_layout: 
function [layout_xy] = get_best_layout(this_dist_tree, layout_seed, options_struct)
	% do multiple HK layouts
	n_layouts 		= 16;
	layout_cell 	= cell(n_layouts, 1);
	energy 			= zeros(n_layouts, 1);
	energy_hk 		= zeros(n_layouts, 1);

	if options_struct.pool_flag
		% set up reproducible rng
		spmd
			cmrg 			= RandStream('mrg32k3a', 'seed', layout_seed);
			RandStream.setGlobalStream(cmrg);
		end

		parfor ii = 1:n_layouts
			% set up random stream
			s 				= RandStream.getGlobalStream();
			s.Substream 	= ii;

			% do layout
			[this_layout, this_energy, this_energy_hk] = harel_koren_layout_faster(this_dist_tree);

			layout_cell{ii}	= this_layout;
			energy(ii)		= this_energy;
			energy_hk(ii)	= this_energy_hk;
		end
	else
		for ii = 1:n_layouts
			% set up random stream
			rng(ii);

			% do layout
			[this_layout, this_energy, this_energy_hk] = harel_koren_layout_faster(this_dist_tree);

			layout_cell{ii}	= this_layout;
			energy(ii)		= this_energy;
			energy_hk(ii)	= this_energy_hk;
		end
	end

	% take best
	[~, best_idx] 	= min(energy);
	layout_xy 		= layout_cell{best_idx};

	% scale layout nicely
	min_xy 			= min(layout_xy, [], 1);
	max_xy 			= max(layout_xy, [], 1);
	layout_xy 		= bsxfun(@rdivide, bsxfun(@minus, layout_xy, min_xy), max_xy - min_xy)*0.8 + 0.1;
end

