%% get_null_distribution: 
function null_distribution = get_null_distribution(n_ref_cells, n_points, n_dims)
	if n_points <= 1000
		fprintf('too few observations for proper null distribution; not significant\n');
		null_distribution 	= Inf;
		return
	end

	% load all null distributions
	load('treetop null distns.mat', 'lookup_table', 'null_cell')

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

	% get this null distribution as outputs
	this_idx 		= find( ...
		lookup_table.n_ref_cells == n_ref_cells_match & ...
		lookup_table.n_points== n_points_match & ...
		lookup_table.n_dims == n_dims_match ...
		);
	if numel(this_idx) ~= 1
		error('null distribution matching went wrong')
	end
	null_distribution 	= null_cell{this_idx};
end
