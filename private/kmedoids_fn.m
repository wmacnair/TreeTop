%% kmedoids_fn: calculates kmedoids algorithm with given distance matrix, and weights
% inputs are: D, distance matrix; kk, number of clusters; weights, vector of weights
% uses algorithm from Park and Jun, 2009
function [cluster_labels, medoids_idx, energy] = kmedoids_fn(D, kk, weights)
	% check inputs
	[weights, nn] 	= check_inputs(D, weights);

	% calculate weighted distance matrix
	weighted_D 		= bsxfun(@times, D, weights);

	% pick starting medoids
	start_idx 			= randsample(nn, kk);

	% initialize
	medoids_idx 		= start_idx;
	[~, cluster_new] 	= min(D(medoids_idx, :), [], 1);

	continue_flag 		= true;

	% loop
	while continue_flag
		% calculate total weighted distance to these medoids
		total_weighted_D 		= weighted_D' * sparse(1:nn, cluster_new, 1, nn, kk, nn);

		% which is best medoid for each of these clusters?
		[~, medoids_idx] 		= min(total_weighted_D,[],1);

		% remember old cluster
		cluster_old 			= cluster_new;

		% assignment step
		[energy, cluster_new] 	= min(D(medoids_idx, :), [], 1);

		% check whether to continue
		continue_flag 			= any(cluster_new ~= cluster_old);
	end
	% tidy up for outputs
	cluster_labels 		= cluster_new;
end

%% check_inputs: 
function [weights, nn] = check_inputs(D, weights)
	nn 	= size(D, 1);

	% are dist matrix and weights appropriate sizes?
	if nn ~= numel(weights)
		error('D and weights must be equal dimensions')
	end

	% is weights a column vector?
	if size(weights, 2) ~= 1
		weights 	= weights';
	end
end

%% calc_v_i: v_i used to find which nodes are closest to the centre of the data
% implemented in a non-memory hungry way
function [v_i] = calc_v_i(D, weighted_D);
	% sum_D 		= sum(D, 2);
	% prop_D 		= bsxfun(@rdivide, weighted_D, sum_D);
	% v_i 		= sum(prop_D, 1);
	
	n_nodes 	= size(D, 1);
	v_i 		= NaN(1, n_nodes);

	sum_D 		= sum(D, 2);
	
	for ii = 1:n_nodes
		prop_D_ii 		= weighted_D(:, ii)./sum_D;
		v_i(ii) 		= sum(prop_D_ii, 1);
	end
end
