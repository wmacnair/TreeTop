%% treetop_pre_run: Run this before running TreeTop, to check that the markers used are useful. 
% Outputs are plots of marginal distributions of all markers, split by input file, and plots of 
% mutual information (MI) between markers. High MI between two markers indicates that they share 
% information, and therefore might be involved in the same process. For each marker, the maximum 
% MI with all other markers is shown; markers with low maximum MI share little information with 
% any other marker, and can be considered for exclusion to improve signal in the data.
function treetop_pre_run(input_struct, options_struct)
	fprintf('\nrunning pre-run analysis for TreeTop\n')

	% get data
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);
	fprintf('1/6 Getting data\n')
	all_struct 		= get_all_files(input_struct);

	% plot marginals
	plot_marginals(all_struct, input_struct, options_struct)

	% calculate and plot MI
	plot_mi(all_struct, input_struct, options_struct)

	fprintf('\ndone.\n\n')
end

%% plot_marginals: 
function plot_marginals(all_struct, input_struct, options_struct)
	% unpack
	all_labels 		= all_struct.all_labels;
	used_data 		= all_struct.used_data;
	used_markers 	= all_struct.used_markers;
	extra_data 		= all_struct.extra_data;
	extra_markers 	= all_struct.extra_markers;

	% identify labels
	unique_labels 	= unique(all_labels);
	label_idx 		= cellfun(@(this_label) strcmp(this_label, all_labels), unique_labels, 'unif', false);

	% plot used_markers
	fprintf('2/6 Plotting marginals of used markers');
	[fig, fig_size] = plot_marginals_once(used_data, used_markers, unique_labels, label_idx, input_struct.used_cofactor);
	plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s used marginals', input_struct.save_stem));
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);

	if ~isempty(extra_data)
		fprintf(', and of extra markers');
		% plot used_markers
		[fig, fig_size] = plot_marginals_once(extra_data, extra_markers, unique_labels, label_idx, input_struct.extra_cofactor);
		plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s extra marginals', input_struct.save_stem));
		plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);
	end
	fprintf('\n')
end

%% plot_marginals_once: 
function [fig, fig_size] = plot_marginals_once(this_data, this_markers, unique_labels, label_idx, cofactor)
	% set up plot
	fig 			= figure('visible', 'off');
	n_plots 		= size(this_data, 2);
	plot_ratio 		= 1.5;
	n_rows 			= ceil(sqrt(n_plots / plot_ratio));
	n_cols 			= ceil(n_plots / n_rows);
	plot_unit 		= 4;
	fig_size 		= [plot_unit*n_cols plot_unit*n_rows];

	if cofactor == 5
		cytof_edges 	= 0:0.5:8;
	else
		cytof_edges 	= [];
	end
	
	% do plots
	for ii = 1:n_plots
		subplot(n_rows, n_cols, ii)
		
		this_col 	= this_data(:, ii);

		hold on
		if isempty(cytof_edges)
			[~, bin_edges] 	= histcounts(this_col, 20);
		else
			bin_edges 		= cytof_edges;
		end
		for jj = 1:size(label_idx, 1)
			histogram(this_col(label_idx{jj}), bin_edges);
		end
		hold off
		xlim([min(bin_edges), max(bin_edges)])
		xlabel('Marker value')
		ylabel('# cells')
		title(this_markers{ii}, 'interpreter', 'none')

		if ii == n_plots
			legend(unique_labels{:})
		end
	end
end

%% plot_mi: plot matrix of MI values, and max / mean MI values by marker
function plot_mi(all_struct, input_struct, options_struct)
	% calculate all MI pairs
	try
		mi_mat 				= calc_mi_mat(all_struct);
	catch lasterr
		fprintf('Seems like the MIToolboxMex is not working. Try compiling it?\nMutual information not calculated.\n');
		return
	end

	% order by max
	max_mi 				= max(mi_mat);
	mean_mi 			= mean(mi_mat);
	n_used 				= size(all_struct.used_data, 2);
	n_total				= size(mi_mat, 1);
	used_idx 			= 1:n_total <= n_used;
	[~, order_idx] 		= sortrows(-[used_idx; max_mi]');

	% tidy up marker names
	all_markers 		= {all_struct.used_markers{:}, all_struct.extra_markers{:}};
	all_markers 		= cellfun(@(str) regexprep(str, '^[0-9]{3}[A-Z][a-z]_', ''), all_markers, 'unif', false);
	all_markers 		= cellfun(@(str) regexprep(str, '_', ' '), all_markers, 'unif', false);

	% put things in right order
	all_markers 		= all_markers(order_idx);
	mi_mat 				= mi_mat(order_idx, order_idx);
	max_mi 				= max_mi(order_idx);
	mean_mi 			= mean_mi(order_idx);

	% plot all MI pairs
	plot_mi_mat(mi_mat, all_markers, n_used, n_total, input_struct, options_struct)

	% plot max and mean values of MI for each marker
	plot_max_mi(mean_mi, max_mi, n_used, n_total, all_markers, input_struct, options_struct)

	% print MI outputs into console
	print_mi_values(n_used, n_total, max_mi, mean_mi, all_markers)
end

%% calc_mi_mat: calculate matrix of MI values between markers
function [mi_mat] = calc_mi_mat(all_struct)

	% calculate MI
	all_data 			= [all_struct.used_data, all_struct.extra_data];
	[n_points, n_total] = size(all_data);
	n_bins 				= ceil(n_points^(1/3));

	% do discretization
	fprintf('3/6 Discretizing data\n')
	discrete_type 		= 'equalwidth';
	switch discrete_type
		case 'equalwidth'
			% find bin edges for each column
			xmin 		= min(all_data, [], 1);
			xmax 		= max(all_data, [], 1);
			xrange 		= xmax - xmin;
			edges 		= arrayfun(@(jj) binpicker(xmin(jj), xmax(jj), n_bins, xrange(jj)/n_bins), 1:n_total, 'unif', false);

			% do discretization
			data_discrete_cell 	= arrayfun(@(jj) discretize(all_data(:, jj), edges{jj}), 1:n_total, 'unif', false);
			data_discrete 		= cell2mat(data_discrete_cell);

		case 'equalfreq'
			% calculate quantiles
			quants 				= arrayfun(@(jj) [-Inf, unique(quantile(all_data(:, jj), n_bins-1)), Inf], 1:n_total, 'unif', false);

			% do discretization
			data_discrete_cell 	= arrayfun(@(jj) discretize(all_data(:, jj), quants{jj}), 1:n_total, 'unif', false);
			data_discrete 		= cell2mat(data_discrete_cell);

		otherwise
			error('invalid discretization type')
	end

	% calculate MI
	fprintf('4/6 Calculating MI\n')
	[mesh_ii, mesh_jj] 	= meshgrid(1:n_total, 1:n_total);
	keep_idx 			= mesh_ii < mesh_jj;
	mesh_ii 			= mesh_ii(keep_idx);
	mesh_jj 			= mesh_jj(keep_idx);

	% plot MI
	mi_mat_long 		= arrayfun(@(ii, jj) mi(data_discrete(:, ii), data_discrete(:, jj)), mesh_ii(:), mesh_jj(:));
	mi_mat 				= zeros(n_total);
	utri_idx 			= sub2ind([n_total, n_total], mesh_ii, mesh_jj);
	mi_mat(utri_idx) 	= mi_mat_long;
	mi_mat 				= mi_mat + mi_mat';
end

%% plot_mi_mat:
function plot_mi_mat(mi_mat, all_markers, n_used, n_total, input_struct, options_struct)
	% do plotting
	fprintf('5/6 Plotting matrix of MI values\n')
	fig 				= figure('visible', 'off');
	imagesc(mi_mat)
	hold on
	ylim_vals 			= ylim();
	n_total 				= size(mi_mat, 1);
	set(gca, 'xtick', 1:n_total);
	set(gca, 'ytick', 1:n_total);

	% remove annoying bit of markers
	set(gca, 'yticklabels', all_markers);
	set(gca, 'xticklabels', all_markers);
	set(gca, 'XTickLabelRotation', -45)

	% if necessary, add line between used and extra markers
	if n_used < n_total
		line([n_used, n_used] + 0.5, ylim_vals, 'linestyle', '-', 'color', 'k');
		line(ylim_vals, [n_used, n_used] + 0.5, 'linestyle', '-', 'color', 'k');
	end
	xlim(ylim_vals)
	ylim(ylim_vals)
	hold off

	% add colorbar in a nice place (position = [left bottom width height])
	plot_pos 	= get(gca, 'Position');
	bar_obj 	= colorbar('Position', [plot_pos(1)+plot_pos(3)+0.01, plot_pos(2) + plot_pos(4)/4, 0.02,  plot_pos(4)/2]);
	ylabel(bar_obj, 'MI (bits)')

	% plot
	plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s MI matrix', input_struct.save_stem));
	fig_size 		= [n_total*0.2*1.2, n_total*0.2];
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);
end

%% plot_max_mi: plots max and mean MI values for each marker
function plot_max_mi(mean_mi, max_mi, n_used, n_total, all_markers, input_struct, options_struct)
	fprintf('6/6 Plotting max and mean MI values\n')

	% set up figure
	fig 				= figure('visible', 'off');
	n_rows 				= 2;
	n_cols 				= 1;

	% plot maximum MI observed per marker
	subplot(n_rows, n_cols, 1)
	plot(1:n_total, mean_mi, '.', 'markersize', 20)
	hold on
	set(gca, 'xtick', 1:n_total);
	set(gca, 'xticklabels', all_markers);
	set(gca, 'XTickLabelRotation', -45)
	ylim_vals 			= ylim();
	ylim_vals(1) 		= 0;

	% if necessary, add divider between used and extra
	if n_used < n_total
		line([n_used, n_used] + 0.5, ylim_vals, 'linestyle', '-', 'color', 'k');
	end

	% label
	ylabel('Mean MI (bits)')
	ylim(ylim_vals)
	hold off

	% plot maximum MI observed per marker
	subplot(n_rows, n_cols, 2)
	plot(1:n_total, max_mi, '.', 'markersize', 20)
	hold on
	set(gca, 'xtick', 1:n_total);
	set(gca, 'xticklabels', all_markers);
	set(gca, 'XTickLabelRotation', -45)
	ylim_vals 			= ylim();
	ylim_vals(1) 		= 0;
	% if necessary, add divider between used and extra
	if n_used < n_total
		line([n_used, n_used] + 0.5, ylim_vals, 'linestyle', '-', 'color', 'k');
	end

	% label
	xlabel('Marker (ordered by max value)')
	ylabel('Maximum MI (bits)')
	ylim(ylim_vals)
	hold off

	% plot
	plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s max MI', input_struct.save_stem));
	fig_size 		= [6, 6];
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);
end

%% print_mi_values: 
function print_mi_values(n_used, n_total, max_mi, mean_mi, all_markers)
	% order MI differently
	type_list 		= [repmat({'used'}, n_used, 1); repmat({'extra'}, n_total - n_used, 1)];
	[~, mi_order] 	= sort(-max_mi);

	% print out MI details
	fprintf('\nList of markers ordered by maximum pairwise MI observed with other markers:\n\n')
	fprintf('Type\tMax MI\tMean MI\tMarker name\n')
	for ii = 1:n_total
		this_idx 	= mi_order(ii);
		fprintf('%s\t%.2f\t%.2f\t%s\n', type_list{this_idx}, max_mi(this_idx), mean_mi(this_idx), all_markers{this_idx});
	end
end
