%% treetop_plots: plot marker values on layout, sample arrangement, and bifurcation outputs
% maybe also do ANOVA thing?
function treetop_plots(input_struct, options_struct)
	% parse inputs
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);

	% get outputs we need
	treetop_struct 		= get_treetop_outputs(input_struct);

	% do / get layout
	layout_struct 		= get_layout_struct(input_struct, options_struct);
	
	% plot marker values nicely (both used and extra markers)
	plot_marker_values_on_treetop(treetop_struct, layout_struct, input_struct, options_struct);

	% plot samples over layout
	plot_samples_on_treetop(treetop_struct, layout_struct, input_struct, options_struct);

	% plot density distribution
	plot_density(treetop_struct, layout_struct, input_struct, options_struct);

	% plot scores for all points over graph layout
	plot_branching_scores(treetop_struct, layout_struct, input_struct, options_struct);

	% plot scores for all points over graph layout
	plot_branch_profiles(treetop_struct, input_struct, options_struct);
end

%% plot_marker_values_on_treetop: 
function [] = plot_marker_values_on_treetop(treetop_struct, layout_struct, input_struct, options_struct)
	% unpack
	layout_xy 		= layout_struct.layout_xy;
	layout_graph 	= layout_struct.layout_graph;

	% define data to loop through
	values_cell 	= {treetop_struct.used_values, 	treetop_struct.extra_values};
	markers_cell 	= {input_struct.used_markers, 	input_struct.extra_markers};
	names_cell 		= {'used', 						'extra'};

	for ii = 1:length(values_cell)
		% unpack
		this_values 	= values_cell{ii};
		this_markers 	= markers_cell{ii};
		this_name 		= names_cell{ii};

		if isempty(this_values)
			continue
		end
		fprintf('plotting mean %s marker values\n', this_name);

		% start figure
		fig 			= figure('Visible','off');
		n_plots 		= size(this_values, 2);
		plot_ratio 		= 1.5;
		n_rows 			= ceil(sqrt(n_plots / plot_ratio));
		n_cols 			= ceil(n_plots / n_rows);
		clim 			= [0, 1];

		% scale marker values
		max_vals 		= max(this_values, [], 1);
		min_vals 		= min(this_values, [], 1);
		scaled_vals = bsxfun( ...
			@times, ...
			bsxfun(@minus, this_values, min_vals), ...
			1./(max_vals - min_vals) ...
			);

		% loop through all selected markers
		for ii = 1:n_plots
			% set up plot
			subplot(n_rows, n_cols, ii)

			% order values
			these_vals 	= scaled_vals(:, ii);
			[~, idx] 		= sort(these_vals);

			% plot
			hold on
			gplot(layout_graph, layout_xy, '-k');
			scatter(layout_xy(idx, 1), layout_xy(idx, 2), [], these_vals(idx), 'filled')
			xlim([0, 1])
			ylim([0, 1])
			caxis(clim);

			% labels
			set(gca, 'XTickLabel', '')
			set(gca, 'YTickLabel', '')
			title(this_markers{ii}, 'interpreter', 'none')
			xlabel('TreeTop 1'); ylabel('TreeTop 2'); 

			hold off
		end

		last_plot 	= get(subplot(n_rows, n_cols, n_cols), 'Position');
		% [left bottom width height]
		colorbar('Position', [last_plot(1)+last_plot(3)+0.01  last_plot(2)  0.02  last_plot(4)])
		% remove white space
		set(gca,'LooseInset', get(gca,'TightInset'));

		% save figure
		plot_stem 	= fullfile(input_struct.output_dir, sprintf('%s %s marker values', input_struct.save_stem, this_name));
		plot_unit 	= 4;
		fig_size 	= [plot_unit*n_cols plot_unit*n_rows];
		plot_fig(fig, plot_stem, options_struct.file_ext, fig_size)
	end
end

%% plot_samples_on_treetop: 
function [] = plot_samples_on_treetop(treetop_struct, layout_struct, input_struct, options_struct)
	fprintf('plotting samples\n');

	% unpack
	layout_xy 		= layout_struct.layout_xy;
	layout_graph 	= layout_struct.layout_graph;
	sample_counts 	= treetop_struct.sample_counts;
	sample_names 	= treetop_struct.sample_names;

	% start figure
	fig 			= figure('Visible','off');
	n_plots 		= size(sample_counts, 2);
	plot_ratio 		= 2;
	n_rows 			= ceil(sqrt(n_plots / plot_ratio));
	n_cols 			= ceil(n_plots / n_rows);

	% scale values
	max_size 		= 500;
	scale_factor 	= max_size / max(sample_counts(:));

	% loop through all selected markers
	for ii = 1:n_plots
		% set up plot
		subplot(n_rows, n_cols, ii)

		% order values
		these_vals 		= scale_factor * sample_counts(:, ii);
		non_zeros 		= these_vals > 0;

		% plot
		hold on
		grey_gplot(layout_graph, layout_xy);
		scatter(layout_xy(non_zeros, 1), layout_xy(non_zeros, 2), these_vals(non_zeros), 'filled');
		xlim([0, 1])
		ylim([0, 1])

		% labels
		set(gca, 'XTickLabel', '')
		set(gca, 'YTickLabel', '')
		xlabel('TreeTop 1')
		ylabel('TreeTop 2')
		title(sample_names{ii}, 'interpreter', 'none')

		hold off
	end

	% save figure
	plot_stem 	= fullfile(input_struct.output_dir, sprintf('%s all samples', input_struct.save_stem));
	plot_unit 	= 4;
	fig_size 	= [plot_unit*n_cols plot_unit*n_rows];
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);

	% only plot largest sample plot if we have more than one sample
	if numel(sample_names) == 1
		return
	end

	% prepare weights
	sample_totals 	= sum(sample_counts, 1);
	sample_props 	= bsxfun(@rdivide, sample_counts, sample_totals);
	[~, max_sample] = max(sample_props, [], 2);

	% do plotting
	fig 			= figure('Visible','off');
	% subplot(2, 3, [1:2, 4:5])
	hold on
	gplot(layout_graph, layout_xy, '-k');
	better_gscatter(layout_xy(:, 1), layout_xy(:, 2), max_sample);
	xlim([0, 1])
	ylim([0, 1])
	hold off

	% add legend, adjust position
	h_all			= findobj(gca,'Type','line');
	h_legend 		= legend(h_all(end-1:-1:1), sample_names, 'fontsize', 6, 'location', 'eastoutside');
	% 	pos_legend 		= get(h_legend, 'position');
	% 	pos_legend(1) 	= 0.9;
	% 	pos_legend(2) 	= 0.4;
	% 	set(h_legend, 'position', pos_legend);

	% labels
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	title('Sample with highest proportion at each point')

	% save figure
	plot_stem 	= fullfile(input_struct.output_dir, sprintf('%s largest sample', input_struct.save_stem));
	fig_size 	= [6, 4];
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size);
end

%% plot_density: mainly for diagnostics on density
function [] = plot_density(treetop_struct, layout_struct, input_struct, options_struct)
	fprintf('plotting density\n');

	% unpack
	layout_xy 			= layout_struct.layout_xy;
	layout_graph 		= layout_struct.layout_graph;
	density 			= treetop_struct.density;
	cell_assignments 	= treetop_struct.cell_assignments;
	celltypes			= treetop_struct.celltype_vector;

	% calculate mean density at each point
	[mean_density, labels]	= grpstats(density, cell_assignments, {'mean', 'gname'});
	% remove any zeros
	if strcmp(labels{1}, '0')
		mean_density 		= mean_density(2:end);
	end

	% start figure
	fig 			= figure('Visible','off');
	n_rows 			= 1;
	n_cols 			= 3;
	plot_ii			= 1;

	% plot mean density at each point
	subplot(n_rows, n_cols, plot_ii)
	plot_ii			= plot_ii + 1;
	hold on
	gplot(layout_graph, layout_xy, '-k');
	scatter(layout_xy(:, 1), layout_xy(:, 2), [], mean_density, 'filled')
	xlim([0, 1])
	ylim([0, 1])

	% labels
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	xlabel('TreeTop 1')
	ylabel('TreeTop 2')
	title({'Mean density at each point', sprintf('(sigma = %.1e)', options_struct.sigma)})
	hold off

	% plot histogram of all density values
	subplot(n_rows, n_cols, plot_ii)
	plot_ii			= plot_ii + 1;
	histogram(density);
	xlim_vals 		= xlim;
	xlim([0, xlim_vals(2)]);
	xlabel('Density')
	ylabel('# cells')
	title('Distribution of all density values')

	% set up celltypes
	celltype_list	= unique(celltypes);
	n_labels		= numel(celltype_list);

	% plot density distributions by label
	subplot(n_rows, n_cols, plot_ii)
	plot_ii			= plot_ii + 1;
	if n_labels > 9
		clr			= jet(n_labels);
	elseif n_labels < 3
		clr			= cbrewer('qual', 'Set1', 3);
		clr 		= clr(1:n_labels, :);
	else
		clr			= cbrewer('qual', 'Set1', n_labels);
	end
	
	hold on
	for kk = 1:n_labels
		% restrict to this label
		this_label	= celltype_list(kk);
		this_idx	= celltypes == this_label;
		
		% plot ECDF
		h(kk)		= cdfplot(density(this_idx));
		set(h(kk), 'Color', clr(kk,:))
	end
	hold off
	% add legend, adjust position
	h_legend 		= legend({char(celltype_list)}, 'FontSize', 4);
	pos_legend 		= get(h_legend,'position');
	pos_legend(1) 	= 0.8;
	pos_legend(2) 	= 0.2;
	set(h_legend, 'position', pos_legend);

	% add other labels
	xlabel('Density')
	ylabel('F(x)')
	title('ECDF of density values by label')
	
	% save figure
	plot_stem 	= fullfile(input_struct.output_dir, sprintf('%s density distribution', input_struct.save_stem));
	plot_unit 	= 4;
	fig_size 	= [plot_unit*n_cols*1.1 plot_unit*n_rows];
	plot_fig(fig, plot_stem, options_struct.file_ext, fig_size)
end

%% plot_branching_scores: 
function [] = plot_branching_scores(treetop_struct, layout_struct, input_struct, options_struct)
	fprintf('plotting branching scores\n');

	% unpack
	branch_scores 		= treetop_struct.branch_scores;
	best_branches 		= treetop_struct.best_branches;
	n_ref_cells 		= length(best_branches);
	[n_points, n_dims] 	= size(treetop_struct.used_data);
	file_ext 			= options_struct.file_ext;
	layout_xy 			= layout_struct.layout_xy;
	layout_graph 		= layout_struct.layout_graph;

	% get cell labelling data
	sample_names 		= treetop_struct.sample_names;
	n_samples 			= length(sample_names);

	% set up figure
	fig 			= figure('visible', 'off');
	if n_samples > 1
		n_cols 		= 4;
	else
		n_cols 		= 3;
	end
	n_rows 			= 1;

	% normalize branch scores with reference to scores from non-branching distribution
	non_branching_distn = get_non_branching_distn(n_ref_cells, n_points, n_dims);
	q_cutoff 			= quantile(non_branching_distn, options_struct.p_cutoff);
	normed_scores 		= branch_scores / q_cutoff;

	% plot scores in increasing order
	subplot(n_rows, n_cols, 1)
	[~, score_idx]	= sort(normed_scores);
	hold on
	gplot(layout_graph, layout_xy, '-k');
	scatter(layout_xy(score_idx, 1), layout_xy(score_idx, 2), 30, normed_scores(score_idx), 'filled');
	xlim([0, 1])
	ylim([0, 1])
	hold off

	% label graph
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	xlabel('TreeTop 1')
	ylabel('TreeTop 2')
	max_score		= max(normed_scores);
	title_str 		= {'Relative branching scores', sprintf('(max = %.1f)', max_score)};
	title(title_str)

	% plot distribution of scores relative to defined cutoff
	subplot(n_rows, n_cols, 2)
	hold on
	ecdf(normed_scores)
	ylim_vals 		= ylim(gca);
	line([1, 1], ylim_vals, 'linestyle', '--', 'color', 'k');
	hold off

	% label graph
	n_branch_pts 	= sum(normed_scores > 1);
	xlabel('Relative branching score')
	ylabel('ECDF(score)')
	title_str 		= {'Distribution of scores', sprintf('(%d higher than non-branching distribution)', n_branch_pts)};
	title(title_str)

	% plot branches for highest score
	subplot(n_rows, n_cols, 3)

	% plot graph as background
	hold on
	grey_gplot(layout_graph, layout_xy);

	% identify non-singleton, non branching point branches
	branch_counts 	= tabulate(best_branches);
	branch_list		= branch_counts(:, 1);
	disp_branches 	= branch_list(branch_counts(:, 1) > 0 & branch_counts(:, 2) > 1);
	show_idx 		= ismember(best_branches, disp_branches);
	if max_score > 1
		better_gscatter(layout_xy(show_idx, 1), layout_xy(show_idx, 2), best_branches(show_idx));
	else
		plot(layout_xy(:, 1), layout_xy(:, 2), '.', 'markersize', 10);
	end

	% plot branching point itself
	split_idx 		= best_branches == 0;
	plot(layout_xy(split_idx, 1), layout_xy(split_idx, 2), '.k', 'markersize', 40);
	xlim([0, 1])
	ylim([0, 1])
	hold off

	% labels
	set(gca, 'XTickLabel', '')
	set(gca, 'YTickLabel', '')
	xlabel('TreeTop 1')
	ylabel('TreeTop 2')
	title('Consensus branches')
	h_all			= findobj(gca,'Type','line');
	branch_names 	= arrayfun(@num2str, disp_branches, 'unif', false);
	h_legend 		= legend(h_all(end-1:-1:2), branch_names, 'fontsize', 6);

	% if worthwhile, plot contingency table of branches vs celltypes
	if n_samples > 1
		% unpack
		cell_assignments 	= treetop_struct.cell_assignments;
		celltype_vector 	= treetop_struct.celltype_vector;

		% remove any outliers
		non_outlier_idx 	= cell_assignments ~= 0;
		cell_assignments 	= cell_assignments(non_outlier_idx);
		celltype_vector 	= celltype_vector(non_outlier_idx);

		% assign branches to all original cells
		branches_by_cell 	= best_branches(cell_assignments);

		% remove singleton branch assignments
		cells_to_keep 		= ismember(branches_by_cell, disp_branches);
		branches_by_cell 	= branches_by_cell(cells_to_keep);
		celltype_vector 	= celltype_vector(cells_to_keep);

		% put in order from top left to top right
		[mean_branch_by_celltype, labels]		= grpstats(branches_by_cell, celltype_vector, {'mean', 'gname'});
		[~, sort_idx]							= sort(mean_branch_by_celltype);
		col_order								= labels(sort_idx);
		
		subplot(n_rows, n_cols, 4)
		plot_contingency_table(branches_by_cell, celltype_vector, [], col_order)
	end

	% save outputs
	plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s branching outputs', input_struct.save_stem));
	plot_unit 		= 4;
	fig_size 		= [plot_unit*n_cols*1.1, plot_unit*n_rows];
	plot_fig(fig, plot_stem, file_ext, fig_size)
end

%% grey_gplot: 
function h2 = grey_gplot(layout_graph, layout_xy)
	gplot(layout_graph, layout_xy, '-k');
	h 				= gca;
	h2 				= get(h, 'Children');
	grey_val 		= 0.8;
	set(h2, 'color', [grey_val, grey_val, grey_val]);
end

%% plot_branch_profiles: plot scores for all points over graph layout
function plot_branch_profiles(treetop_struct, input_struct, options_struct)
	if isfield(treetop_struct, 'mean_tree_dist')
		fprintf('plotting markers by branch\n')
	else
		fprintf('MST distance to branch point not calculated; not plotting markers by branch\n')
		return
	end

	% unpack
	mean_tree_dist 	= treetop_struct.mean_tree_dist;
	best_branches 	= treetop_struct.best_branches;

	% set up branch vars
	branch_point 	= find(best_branches == 0);
	unique_branches = setdiff(unique(best_branches), 0);
	n_branches 		= numel(unique_branches);
	palette			= cbrewer('qual', 'Set1', n_branches);

	% get ordering for each branch
	branch_struct_cell 	= arrayfun(@(kk) get_branch_order(best_branches, mean_tree_dist, kk), 1:n_branches, 'unif', false);

	% loop through used, extra
	values_cell 	= {treetop_struct.used_values, treetop_struct.extra_values};
	plot_names 		= {'used', 'extra'};
	markers_cell 	= {input_struct.used_markers, input_struct.extra_markers};
	anova_cell 		= {treetop_struct.anova_used, treetop_struct.anova_extra};

	for ii = 1:2
		% unpack
		sel_values 		= values_cell{ii};
		sel_name 		= plot_names{ii};
		sel_markers 	= markers_cell{ii};
		sel_anova 		= anova_cell{ii};
		if isempty(sel_values)
			fprintf('no %s markers; skipping\n', sel_name)
			continue
		end

		% set up figure
		fig 			= figure('visible', 'off');
		plot_ratio 		= 1.5;
		n_plots			= size(sel_values, 2);
		n_rows 			= ceil(sqrt(n_plots / plot_ratio));
		n_cols 			= ceil(n_plots / n_rows);
		
		% order markers by anova
		[~, anova_idx] 	= sort(sel_anova);

		% separate plot for each marker
		for jj = 1:n_plots
			subplot(n_rows, n_cols, jj)
			hold on

			% do in order of biggest differences
			marker_idx		= anova_idx(jj);

			% get branch_1 values
			% report predictions for branches 1 and kk
			% store all branch_1 predictions
			% define variable for branch_1_predictions
			smooth_1_all 	= zeros(n_branches-1, sum(branch_struct_cell{1}.this_branch));
			int_lo_1_all 	= zeros(n_branches-1, sum(branch_struct_cell{1}.this_branch));
			int_hi_1_all 	= zeros(n_branches-1, sum(branch_struct_cell{1}.this_branch));

			% separate line for each branch
			for kk = 2:n_branches
				% calculate smoothed values along branch
				[smooth_1, int_1, smooth_kk, int_kk] 	= get_branch_vals(branch_struct_cell, kk, sel_values, marker_idx);
				smooth_1_all(kk-1, :) 					= smooth_1;
				int_lo_1_all(kk-1, :) 					= int_1(:, 1);
				int_hi_1_all(kk-1, :) 					= int_1(:, 2);

				% set up colour
				branch_col 				= palette(kk, :);

				% % plot original values
				% plot(sorted_dist, branch_vals, '.', 'color', branch_col);

				% put smoothed values in right order
				sort_dist_kk 			= branch_struct_cell{kk}.sorted_dist;
				plot(sort_dist_kk, smooth_kk, '-', 'color', branch_col, 'linewidth', 2);
				plot(sort_dist_kk, smooth_kk - int_kk(:,1), ':', 'color', branch_col, 'linewidth', 1);
				plot(sort_dist_kk, smooth_kk + int_kk(:,2), ':', 'color', branch_col, 'linewidth', 1);
			end

			% plot branch 1
			smooth_1_mean 		= mean(smooth_1_all, 1);
			int_lo_1_mean 		= mean(int_lo_1_all, 1);
			int_hi_1_mean 		= mean(int_hi_1_all, 1);
			branch_col 			= palette(1, :);
			sort_dist_1 		= branch_struct_cell{1}.sorted_dist;
			plot(sort_dist_1, smooth_1_mean, '-', 'color', branch_col, 'linewidth', 2);
			plot(sort_dist_1, smooth_1_mean - int_lo_1_mean, ':', 'color', branch_col, 'linewidth', 1);
			plot(sort_dist_1, smooth_1_mean + int_hi_1_mean, ':', 'color', branch_col, 'linewidth', 1);

			% % plot cutpoint itself
			% cutpoint_val 	= sel_values(branch_point, marker_idx);
			% plot(0, cutpoint_val, '.', 'color', 'k', 'markersize', 10)

			% tidy up plot
			ylabel('Marker value')
			xlabel('Distance from branching point')
			title_str		= sprintf('%s (p = %.1e)', sel_markers{marker_idx}, sel_anova(marker_idx));
			title(title_str, 'interpreter', 'none')
			hold off

			% ylim([0,1])
		end

		% save outputs
		plot_stem 		= fullfile(input_struct.output_dir, sprintf('%s %s marker profiles', input_struct.save_stem, sel_name));
		plot_unit 		= 4;
		fig_size 		= [plot_unit*n_cols*1.1, plot_unit*n_rows];
		plot_fig(fig, plot_stem, options_struct.file_ext, fig_size)
	end
end

%% get_branch_order: 
function [branch_struct] = get_branch_order(best_branches, mean_tree_dist, jj)
	% extract data
	this_branch 				= best_branches == jj | best_branches == 0;
	branch_dist 				= mean_tree_dist(this_branch);
	[sorted_dist, branch_order] = sort(branch_dist);

	% put into struct
	branch_struct 	= struct( ...
		'this_branch', 	this_branch, ...
		'branch_dist', 	branch_dist, ...
		'branch_order', branch_order, ...
		'sorted_dist', 	sorted_dist ...
		);
end

%% get_branch_vals: 
function [smooth_1, int_1, smooth_kk, int_kk] = get_branch_vals(branch_struct_cell, kk, sel_values, marker_idx)
	% get values along branch 1 in right order
	branch_1_struct 	= branch_struct_cell{1};
	branch_1 			= branch_1_struct.this_branch;
	order_1 			= branch_1_struct.branch_order;
	sort_dist_1 		= branch_1_struct.sorted_dist;
	vals_1 				= sel_values(branch_1, marker_idx);
	vals_1 				= vals_1(order_1);

	% get values along branch kk in right order
	branch_kk_struct 	= branch_struct_cell{kk};
	branch_kk 			= branch_kk_struct.this_branch;
	order_kk 			= branch_kk_struct.branch_order;
	sort_dist_kk 		= branch_kk_struct.sorted_dist;
	vals_kk 			= sel_values(branch_kk, marker_idx);
	vals_kk 			= vals_kk(order_kk);

	% do smoothed fit to both branches together
	% gp_fit 			= fitrgp(sorted_dist(:), branch_vals(:), 'fitmethod', 'fic', 'predictmethod', 'fic');
	fit_obj 			= fitrgp([sort_dist_1(:); sort_dist_kk(:)], [vals_1(:); vals_kk(:)]);

	% do separate predictions
	[smooth_1, int_1] 	= predict(fit_obj, sort_dist_1');
	[smooth_kk, int_kk] = predict(fit_obj, sort_dist_kk');

	% double up
	int_1 				= [int_1, int_1];
	int_kk 				= [int_kk, int_kk];

	% alternative smoothers considered:
	% smooth_vals 	= smooth(sorted_dist, branch_vals, 50, 'rlowess');
	% fit_obj 		= fit([sorted_dist(:); 0], [branch_vals(:); point_val], 'poly2');
	% smooth_vals 	= feval(fit_obj, [sorted_dist(:); 0]);
end
