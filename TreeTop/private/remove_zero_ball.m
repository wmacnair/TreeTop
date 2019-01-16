%% remove_zero_ball: removes 'ball' of observations which are close to zero
function [all_struct] = remove_zero_ball(all_struct, options_struct)
	if ~isfield(options_struct, 'zero_ball_flag') || options_struct.zero_ball_flag == false
		return
	end

	% calculate L1 values
	used_data 		= all_struct.used_data;
	switch options_struct.metric_name
		case {'L1', 'angle'}
			dist_vals 	= sum(abs(used_data), 2);
		case 'L2'
			dist_vals 	= sum(used_data.^2, 2);
		otherwise
			error('haven''t implemented this distance yet...')
	end

	% apply cutoff
	cutoff_val 				= quantile(dist_vals, 0.01);
	keep_idx 				= dist_vals > cutoff_val;
	all_struct.used_data 	= all_struct.used_data(keep_idx, :);
	all_struct.all_labels 	= all_struct.all_labels(keep_idx);
	if ~isempty(all_struct.extra_data)
		all_struct.extra_data 	= all_struct.extra_data(keep_idx, :);
	end
	if isfield(all_struct, 'cell_assignments_top')
		all_struct.cell_assignments_top 	= all_struct.cell_assignments_top(keep_idx, :);
	end
end
