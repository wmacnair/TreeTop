%% get_graph_struct: 
function [graph_struct] = get_graph_struct(input_struct)
	fprintf('loading summaries of ensemble of trees\n')
	
	% open needed graphs
	[inv_freq_graph, dist_graph] 				= get_tree_files(input_struct);

	% calculate maximally sparse connected graphs for each of these
	[sparse_inv_freq_graph, sparse_dist_graph] = calc_sparse_graphs(inv_freq_graph, dist_graph);

	% calculate graphs which are 10% sparser than full graphs
	[q10_inv_freq_graph, q10_dist_graph] 		= calc_10_percent_graphs(inv_freq_graph, dist_graph);

	% put all six into graph structure, with names
	graph_struct = struct( ...
		'inv_adj_matrix', 	{dist_graph, q10_dist_graph, sparse_dist_graph, inv_freq_graph, q10_inv_freq_graph, sparse_inv_freq_graph}, ...
		'name', 			{'Distance graph', 'Semi-sparse distance graph', 'Maximally sparse distance graph', 'Freq graph', 'Semi-sparse freq graph', 'Maximally sparse freq graph'}, ...
		'file_suffix', 		{'_full_dist_layout', '_semi_dist_layout', '_sparse_dist_layout', '_full_freq_layout', '_semi_freq_layout', '_sparse_freq_layout'} ...
		);

	% apply inversion to all graphs, then add outputs back into the structure
	adj_matrix 			= cellfun(@invert_sparse_graph, {graph_struct.inv_adj_matrix}, 'Unif', false);
	[graph_struct(:).adj_matrix] = adj_matrix{:};

	% apply booleanisation to all graphs, then add outputs back into the structure
	boolean_adj_matrix 	= cellfun(@booleanise_sparse_graph, {graph_struct.inv_adj_matrix}, 'Unif', false);
	[graph_struct(:).boolean_adj_matrix] = boolean_adj_matrix{:};
end

%% get_tree_files: 
function [inv_freq_graph, union_graph] = get_tree_files(input_struct)
	% unpack
	output_dir 			= input_struct.output_dir;
	save_stem 			= input_struct.save_stem;

	% open tree based on frequency of edges
	inv_freq_graph_file = fullfile(output_dir, sprintf('%s_freq_union_tree.mat', save_stem));
	inv_freq_graph 		= sparse(importdata(inv_freq_graph_file));

	% open tree based on mean treeSNE distances
	union_graph_file 	= fullfile(output_dir, sprintf('%s_union_tree.mat', save_stem));
	union_graph 		= sparse(importdata(union_graph_file));
end

%% calc_10_percent_graphs: calculates graphs with bottom 10% frequency edges removed
function [q10_inv_freq_graph, q10_union_graph] = calc_10_percent_graphs(inv_freq_graph, union_graph)
	% get info from inverse frequency graph
	[ii jj ss] 	= find(inv_freq_graph);
	[mm nn] 	= size(inv_freq_graph);

	% extract quantiles of frequencies
	[cdf_quantiles cdf_inv_freqs] = ecdf(ss);
	% find which corresponds best to 10% (errs on less sparse side)
	this_quantile 		= 0.9;
	ecdf_quantile_idx 	= max( find(this_quantile > cdf_quantiles) );
			
	% specify quantile and inverse frequency corresponding to 10% less
	ecdf_quantile 		= cdf_quantiles(ecdf_quantile_idx);
	cdf_inv_freq 		= cdf_inv_freqs(ecdf_quantile_idx);

	% mess about to turn this into a new graph
	keep_idx 			= ss <= cdf_inv_freq;
	sparse_graph 		= sparse(ii(keep_idx), jj(keep_idx), ss(keep_idx), mm, nn);

	% how many components in this test graph?
	no_components 		= max(components(sparse_graph));
	% count how many components, only true if there is only one components
	is_connected 		= no_components == 1;

	if ~is_connected
		q10_inv_freq_graph 	= [];
		q10_union_graph 	= [];
	else
		% disp(['10 percent graph corresponds to cutoff of ' num2str(cdf_inv_freq) ' and exclusion of ' num2str(round((1-ecdf_quantile)*100)) '% of the edges']);

		[union_ii union_jj union_ss] 	= find(union_graph);
		% check that graphs are same as above
		if ~and(isequal(ii, union_ii), isequal(jj, union_jj))
			error('union graph and inv freq graph don''t have matching non-zero locations')
		end

		% put sparse union graph together
		sparse_union_graph 	= sparse(union_ii(keep_idx), union_jj(keep_idx), union_ss(keep_idx), mm, nn);

		% make outputs
		q10_inv_freq_graph 	= sparse_graph;
		q10_union_graph 	= sparse_union_graph;
	end
end

%% calc_sparse_graphs: calculates maximally sparse graphs from the inputs, by removing edges with low
% frequency until it is not possible to exclude more without disconnecting the graph
function [sparse_inv_freq_graph, sparse_union_graph] = calc_sparse_graphs(inv_freq_graph, union_graph)
	% test cases:
	% - graph is unconnected: 
	% nn = 100; A = magic(nn); A = A + A'; A(:,1) = 0; A(1,:) = 0; A(1:nn+1:nn*nn) = 0;
	% inv_freq_graph = sparse(A);
	% - graph is connected even at maximal freq requirement
	% nn = 100; A = magic(nn); A = A + A'; A(:,1) = 1; A(1,:) = 1; A(1:nn+1:nn*nn) = 0;
	% inv_freq_graph = sparse(A);

	% get list of inverse frequencies
	inv_freq_list = full(unique(inv_freq_graph));

	% set initial unconnected_freq and connected_freq
	unconnected_freq_idx = 1; connected_freq_idx = numel(inv_freq_list) + 1;

	% have we found the cutoff yet?
	while connected_freq_idx - unconnected_freq_idx > 1
		% set test_freq
		test_freq_idx = get_test_idx(unconnected_freq_idx, connected_freq_idx);

		% test for connectedness
		% is_connected = test_connectedness(inv_freq_graph, inv_freq_list, test_freq_idx);
		is_connected = test_proposed_idx(inv_freq_graph,  inv_freq_list, test_freq_idx);

		% update connected / unconnected frequencies
		if is_connected
			connected_freq_idx = test_freq_idx;
		else
			unconnected_freq_idx = test_freq_idx;
		end
	end

	% use identified connected_freq_idx to define output graphs
	if connected_freq_idx > numel(inv_freq_list)
		error('Graph not connected');
	else
		freq_boolean 			= inv_freq_graph <= inv_freq_list(connected_freq_idx);
		sparse_inv_freq_graph 	= inv_freq_graph .* freq_boolean;
		sparse_union_graph 		= union_graph .* freq_boolean;
	end
end

%% invert_sparse_graph: takes sparse graph as input, returns sparse graph with inverted values for non-zero entries
function [inverted_graph] = invert_sparse_graph(input_graph)
	if ~issparse(input_graph)
		error('input_graph must be sparse')
	end
	% get original graph entries
	[ii jj ss] 	= find(input_graph);
	[mm nn] 	= size(input_graph);
	% invert nonzeros
	inverted_ss = 1./ss;
	% make new graph
	inverted_graph = sparse(ii, jj, inverted_ss, mm, nn);
end

%% booleanise_sparse_graph: takes sparse graph as input, returns sparse graph with inverted values for non-zero entries
function [boolean_graph] = booleanise_sparse_graph(input_graph)
	if ~issparse(input_graph)
		error('input_graph must be sparse')
	end
	% get original graph entries
	[ii jj ss] 		= find(input_graph);
	[mm nn] 		= size(input_graph);
	% invert nonzeros
	boolean_ss 		= repmat(1, size(ss));
	% make new graph
	boolean_graph 	= sparse(ii, jj, boolean_ss, mm, nn);
end

%% get_test_freq_idx: finds midpoint between two input indices
function [test_idx] = get_test_idx(unconnected_idx, connected_idx)
	test_idx 		= floor((connected_idx - unconnected_idx)/2) + unconnected_idx;
end

%% test_connectedness: checks whether a graph is connected for a given test frequency
function [is_connected] = test_connectedness(inv_freq_graph, inv_freq_list, test_freq_idx)
	% define graph to test by setting to zero values with too high inverse freqencies
	[ii jj ss] 		= find(inv_freq_graph);
	[mm nn] 		= size(inv_freq_graph);
	sparse_ss 		= ss <= inv_freq_list(test_freq_idx);
	test_graph 		= sparse(ii, jj, sparse_ss, mm, nn);

	% how many components in this test graph?
	no_components 	= max(components(test_graph));
	% count how many components, only true if there is only one components
	is_connected 	= no_components == 1;
end

%% test_proposed_idx: checks whether a graph is connected for a given test frequency
function [is_connected] = test_proposed_idx(input_graph, edge_value_list, test_edge_value_idx, varargin)
	% define which direction to test in
	if nargin == 3
		test_le = true;
	elseif nargin == 4
		test_le = varargin{1};
	else
		error('arg')
	end

	% define graph to test by setting to zero values with too high inverse freqencies
	[ii jj ss] 		= find(input_graph);
	[mm nn] 		= size(input_graph);
	if test_le
		sparse_ss 	= ss <= edge_value_list(test_edge_value_idx);
	else
		sparse_ss 	= ss >= edge_value_list(test_edge_value_idx);
	end
	test_graph 		= sparse(ii, jj, sparse_ss, mm, nn);

	% check whether this graph is connected
	is_connected 	= test_connected(test_graph);
end

%% test_connected: 
function [is_connected] = test_connected(test_graph)
	% how many components in this test graph?
	no_components = max(components(test_graph));
	% count how many components, only true if there is only one components
	is_connected = no_components == 1;
end

% %% report_sparsity: 
% function [] = report_sparsity(graph_struct)
% 	% get graphs
% 	graphs_cell 	= {graph_struct.inv_adj_matrix};
% 	names_cell 		= {graph_struct.name};

% 	% calculate sparsity
% 	non_zeros 		= cellfun(@nnz, graphs_cell);
% 	total_size 		= cellfun(@numel, graphs_cell);
% 	sparsity 		= non_zeros ./ total_size;

% 	% display
% 	max_name_length = max(cellfun(@length, names_cell));
% 	name_col_length = max_name_length + 5;

% 	fprintf('sparsity of outputs:\n');
% 	spacer 			= horzcat(repmat(' ', 1, max_name_length - length('graph')));
% 	fprintf('graph%snnz\tsize\tsparsity\n', spacer);

% 	for ii = 1:numel(graph_struct)
% 		this_name 	= names_cell{ii};
% 		name_length = length(this_name);
% 		spacer 		= horzcat(repmat(' ', 1, max_name_length - name_length));
% 		fprintf('%s%s%d\t%d\t%.2f%%\n', this_name, spacer, non_zeros(ii), total_size(ii), sparsity(ii)*100);
% 	end
% end
