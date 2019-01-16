%% plot_contingency_table: 
function plot_contingency_table(branches, celltypes, row_order, col_order)
	% do crosstab of branches against gates
	[branch_xtab, ~, ~, labels]		= crosstab(branches, celltypes);

	% calculate NMI between them
	[~, ~, branch_int] 	= unique(branches);
	[~, ~, cell_int] 	= unique(celltypes);
	this_nmi			= nmi(branch_int, cell_int);

	% mess about with labels
	row_labels 		= labels(:, 1);
	row_labels 		= row_labels(~cellfun(@isempty, row_labels));
	col_labels 		= labels(:, 2);
	col_labels 		= col_labels(~cellfun(@isempty, col_labels));

	% normalize by celltype
	normed_xtab 	= bsxfun(@rdivide, branch_xtab, sum(branch_xtab, 1));
	other_normed 	= bsxfun(@rdivide, branch_xtab, sum(branch_xtab, 2));

	% % % put clusters in right order
	% % [~, row_order] 	= sort(row_labels);
	% % row_labels 		= row_labels(row_order);
	% % normed_xtab 	= normed_xtab(row_order, :);

	% % get right order for celltypes in right order
	% link 			= linkage(normed_xtab', 'complete');
	% dist_mat 		= pdist(normed_xtab');
	% col_order 		= optimalleaforder(link, dist_mat);

	% % get right order for clusters
	% link 			= linkage(other_normed, 'complete');
	% dist_mat 		= pdist(other_normed);
	% row_order 		= optimalleaforder(link, dist_mat);

	% % put everything in order
	% row_labels 		= row_labels(row_order);
	% col_labels 		= col_labels(col_order);
	% normed_xtab 	= normed_xtab(row_order, col_order);
	if exist('row_order', 'var') && ~isempty(row_order)
		if ~isequal(sort(row_labels(:)), sort(row_order(:)))
			error('row_order does not match row labels')
		end
		row_idx 		= cellfun(@(this_label) find(strcmp(this_label, row_labels)), row_order);
		row_labels 		= row_labels(row_idx);
		normed_xtab 	= normed_xtab(row_idx, :);
	end
	if exist('col_order', 'var') && ~isempty(col_order)
		if ~isequal(sort(col_labels(:)), sort(col_order(:)))
			error('col_order does not match col labels')
		end
		col_idx 		= cellfun(@(this_label) find(strcmp(this_label, col_labels)), col_order);
		col_labels 		= col_labels(col_idx);
		normed_xtab 	= normed_xtab(:, col_idx);
	end

	% add text of NMI over top
	val_labels 			= arrayfun(@(x) sprintf('%.1f', x), normed_xtab, 'unif', false);
	[mesh_x, mesh_y] 	= meshgrid(1:size(normed_xtab, 2), 1:size(normed_xtab, 1));
	non_zero_idx 		= normed_xtab(:) > 0;

	% plot heatmap of celltypes per cluster, normalized by celltype total (?)
	imagesc(normed_xtab)
	text(mesh_x(non_zero_idx), mesh_y(non_zero_idx), val_labels(non_zero_idx), 'FontSize', 8, 'HorizontalAlignment', 'center')

	% add labels
	ax 				= gca;
	set(gca,'TickLabelInterpreter','none')
	col_ticks 		= 1:size(normed_xtab, 2);
	xlim([min(col_ticks)-0.5, max(col_ticks)+0.5])
	set(ax, 'XTick', col_ticks);
	set(ax, 'XTickLabel', col_labels);
	ax.XTickLabelRotation 	= 45;
	row_ticks 		= 1:size(normed_xtab, 1);
	set(ax, 'YTick', row_ticks);
	set(ax, 'YTickLabel', row_labels);

	% add title
	title({'Allocation of gates to branches', sprintf('(NMI = %.2f)', this_nmi)})
end
