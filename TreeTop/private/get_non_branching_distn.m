%% get_non_branching_distn: 
function non_branching_distn = get_non_branching_distn(n_ref_cells, n_points, n_dims)
	if n_points <= 1000
		fprintf('too few observations for proper non-branching comparison distribution; all scores normalized to 0\n');
		non_branching_distn 	= Inf;
		return
	end

	% load all non-branching distributions
	load('treetop non-branching distns.mat', 'lookup_table', 'non_branching_cell')

	% find closest distribution
	n_ref_cells_list 	= unique(lookup_table.n_ref_cells);
	n_points_list 		= unique(lookup_table.n_points);
	n_dims_list 		= unique(lookup_table.n_dims);

	% find which are closest
	ref_cells_idx = max(find(n_ref_cells_list <= n_ref_cells));
	if isempty(ref_cells_idx)
		ref_cells_idx = 1;
	end
	n_ref_cells_match	= n_ref_cells_list(ref_cells_idx);
	points_idx 		= max(find(n_points_list <= n_points));
	if isempty(points_idx)
		points_idx 		= 1;
	end
	n_points_match	= n_points_list(points_idx);
	dims_idx 		= max(find(n_dims_list <= n_dims));
	if isempty(dims_idx)
		dims_idx 		= 1;
	end
	n_dims_match	= n_dims_list(dims_idx);

	% get this non-branching distribution as outputs
	this_idx 		= find( ...
		lookup_table.n_ref_cells == n_ref_cells_match & ...
		lookup_table.n_points== n_points_match & ...
		lookup_table.n_dims == n_dims_match ...
		);
	if numel(this_idx) ~= 1
		error('non-branching distribution matching went wrong')
	end
	non_branching_distn 	= non_branching_cell{this_idx};
end
