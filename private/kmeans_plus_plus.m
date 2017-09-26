% kmeans_plus_plus.m
% See Scalable K-Means++, Bahmani et al., 2012
% Initializes first cluster selection for k-means, using squared distances from other points to ensure points are well spaced 
% Approximate but distributed method.
% X		rows = observations, columns = fields
% kk 	number of clusters
% ll 	oversampling factor (default is to equal k)
% rr 	number of rounds to run
function [cluster_labels, centroid_idx] = kmeans_plus_plus(X, kk, ll, rr, options, paramset)
	% check inputs ok
	parameters 		= check_inputs(X, kk, ll, rr, options);

	% partition X as evenly as possible amongst cores
	partition_idx 	= calc_partition_idx(parameters);
	X_split 		= mat2cell(X, partition_idx, parameters.n_var);

	continue_flag 	= true;
	while continue_flag
		% do clustering, but check that right number of entries is coming out
		[cluster_labels, centroid_idx] 	= do_clustering(parameters, X, X_split, partition_idx, paramset);

		% check that we have right number, if not, repeat
		n_unique_centroids 				= numel(unique(centroid_idx));
		if n_unique_centroids == kk
			continue_flag 	= false;
		else
			fprintf('didn''t find kk points; repeating seeding\n')
		end
	end
end

%% check_inputs: 
function [parameters] = check_inputs(X, kk, ll, rr, options)
	% verbosity flag
	if ~isfield(options, 'verbose')
		verbose 	= false;
	else
		verbose 	= options.verbose;
	end

	% check what distance measure to use
	if ~isfield(options, 'metric')
		metric = 'L1';
	else
		metric = options.metric;
	end

	% check whether we should use parfor or for
	if ~isfield(options, 'pool_flag')
		% default to no
		pool_flag 	= false;
	else
		pool_flag 	= options.pool_flag;
	end

	% set default number of cores
	n_cores 	= 1;

	% if we're running with a pool, check that the pool actually exists
	if pool_flag
		matlab_version 	= version;

		switch matlab_version
			case '8.1.0.604 (R2013a)'
				pool_size 	= matlabpool('size');
				if pool_size == 0
					warning('No pool; run done without parallelisation.')
					pool_flag 	= false;
				else
					n_cores 	= pool_size;
				end

			otherwise
				current_pool 	= gcp('nocreate'); % If no pool, do not create new one.
				% n_cores = matlabpool('size');
				if isempty(current_pool)
					warning('no pool present; running in serial')
					pool_flag 	= false;
				else
					% if pool does exist, use this as number of cores
					n_cores 	= current_pool.NumWorkers;
				end
		end
	end

	% check whether we have enough observations
	if size(X, 1) < kk
		error('Fewer observations than clusters');
	end

	% 
	% ensure that we get enough values in each iteration
	if ll * rr < kk
		rr 	= max(ceil(kk / ll), 5);
		warning(['Insufficient rounds, value of ', int2str(rr), ' used instead (= ceiling(k/l) )']);
	end

	% set up variables
	[n_obs, n_var] 		= size(X);
	parameters 			= struct();

	parameters.kk 			= kk;
	parameters.rr 			= rr;
	parameters.ll 			= ll;
	parameters.n_cores 		= n_cores;
	parameters.n_obs 		= n_obs;
	parameters.n_var 		= n_var;
	parameters.metric 		= metric;
	parameters.pool_flag 	= pool_flag;
	parameters.verbose 		= verbose;
end

%% partition_idx: 
function [partition_idx] = calc_partition_idx(parameters)
	% define useful variables
	n_obs 			= parameters.n_obs;
	n_cores 		= parameters.n_cores;

	% calculate requirements for partition
	n_per_partition = floor(n_obs / n_cores);
	n_remainder 	= mod(n_obs, n_cores);

	% make partition
	partition_idx 					= ones(n_cores,1) * n_per_partition;
	partition_idx(1:n_remainder) 	= partition_idx(1:n_remainder) + 1;
end

%% do_clustering: have this separate so can repeat if necessary
function [cluster_labels, centroids_idx] = do_clustering(parameters, X, X_split, partition_idx, paramset)
	% unpack
	kk 				= parameters.kk;
	rr 				= parameters.rr;

	% pick first point uniformly at random
	raw_centroids_idx 	= randsample(parameters.n_obs, 1);
	% use this to define centroid_X
	centroid_X 			= X(raw_centroids_idx, :);

	% keep adding new points until we have done at least rr rounds, and we have at least k points
	% check that # centroid_X is at least k
	fprintf('k-means++ to identify reference nodes');
	continue_flag 	= true;
	ii 				= 1;
	while continue_flag
		if ii <= rr
			fprintf('\nextra round to ensure sufficient nodes');
		end

		% calculate induced probability distribution (in distributed way)
		phi 			= par_calc_distn(X_split, centroid_X, parameters);

		% pick next points to be included and update raw_centroids_idx
		[raw_centroids_idx, centroid_X] 	= next_points(X, raw_centroids_idx, phi);

		% loop admin
		ii 				= ii + 1;
		n_centroids 	= numel(raw_centroids_idx);

		% want to carry on until we have reached the right number of rounds, and have sufficient centroid_X
		continue_flag 	= ii <= rr | n_centroids < kk;
	end
	fprintf('\n');

	% recluster points into k clusters
	centroids_idx 		= cluster_into_k(X, X_split, raw_centroids_idx, parameters, paramset);

	% pick centroid closest to cluster_centroids
	cluster_labels 		= get_cluster_labels(X, X_split, centroids_idx, parameters);
end

%% par_calc_distn: calculates distribution in distributed way
function [phi] = par_calc_distn(X_split, centroid_X, parameters)
	% define useful variable
	n_cores 	= parameters.n_cores;

	% define storage variable
	phi_split 	= cell(n_cores, 1);
	
	% calculate distribution for each individual section
	if parameters.pool_flag
		parfor ii = 1:n_cores
			phi_split{ii} 	= calc_distn(X_split{ii}, centroid_X, parameters);
		end
	else
		for ii = 1:n_cores
			phi_split{ii} 	= calc_distn(X_split{ii}, centroid_X, parameters);
		end
	end

	% combine split phis back into one, and add upscaling factor
	phi 		= cell2mat(phi_split);
	phi 		= parameters.ll * phi / sum(phi);
end

%% calc_distn: 
function [phi_split] = calc_distn(X, centroid_X, parameters)
	if parameters.verbose
		show_var_size(size(X,1), 1);
	end

	switch parameters.metric
		case 'L1'
			% calculate distances between points in X and points in the set of centroid_X ()
			phi_split 	= comp_phi_L1(X, centroid_X);
		case 'L2'
			% calculate distances between points in X and points in the set of centroid_X ()
			phi_split 	= comp_phi_euclidean(X, centroid_X);
		case 'angle'
			% calculate distances between points in X and points in the set of centroid_X ()
			phi_split 	= comp_phi_angle(X, centroid_X);
		case 'corr'
			% calculate distances between points in X and points in the set of centroid_X ()
			phi_split 	= comp_phi_corr(X, centroid_X);
		case 'abs_corr'
			% calculate distances between points in X and points in the set of centroid_X ()
			phi_split 	= comp_phi_abs_corr(X, centroid_X);
		otherwise
			error('metric not recognised')
	end
end

%% next_points: selects next points in non-distributed way
function [C, centroid_X] = next_points(X, C, phi)
	% do random sample of uniform variable
	rand_sample 	= rand(numel(phi), 1);
	
	% compare to distribution to decide which were selected
	cluster_idx 	= find(rand_sample < phi);

	% add new cluster points
	C 			= union(C, cluster_idx);
	centroid_X 	= X(C, :);
end

%% cluster_into_k: 
function [centroids_idx] 	= cluster_into_k(X, X_split, raw_centroids_idx, parameters, paramset)
	fprintf('choosing k centroids from initial points\n')
	% not_done_yet 	= true;
	% while not_done_yet

	% unpack
	kk 					= parameters.kk;

	% where are the centroids?
	centroid_X 			= X(raw_centroids_idx, :);
	% calculate weights for each centroids
	% centroid_weights 	= par_calc_centroid_weights(X_split, centroid_X, parameters);
	centroid_weights 	= ones(size(raw_centroids_idx));

	if parameters.verbose
		show_var_size(size(centroid_X,1),size(centroid_X,1));
	end

	% remove some unnecessary variables to improve memory
	clear X X_split

	% calc appropriate distance matrix
	switch parameters.metric
		case 'L1'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_L1(centroid_X, centroid_X);
		case 'L2'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_euclidean(centroid_X, centroid_X);
		case 'angle'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_angle(centroid_X, centroid_X);
		case 'corr'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_corr(centroid_X, centroid_X);
		case 'abs_corr'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_abs_corr(centroid_X, centroid_X);
		otherwise
			error('metric not recognised')
	end

	% remove some unnecessary variables to improve memory
	clear centroid_X

	% calculate weighted cluster
	[~, centroids_idx_idx] 	= kmedoids_fn(D, kk, centroid_weights);

	% which of the original centroids does this correspond to?
	centroids_idx 			= raw_centroids_idx(centroids_idx_idx);

	% have we got the right number of clusters?
	if numel(centroids_idx) ~= kk
		error('got wrong number of clusters')
	end
end

%% par_calc_centroid_weights: calculates how many cells each centroid represents, in a distributed way
function [centroid_weights] = par_calc_centroid_weights(X_split, centroid_X, parameters)
	% define useful variable
	n_cores = parameters.n_cores;

	% define storage variable
	centroid_weights_split = cell(n_cores, 1);

	% choose parallel or not
	if parameters.pool_flag
		% find how many points in X are closest to each centroid
		parfor ii = 1:n_cores
			centroid_weights_split{ii} 	= calc_centroid_weights(X_split{ii}, centroid_X, parameters);
		end
	else
		% find how many points in X are closest to each centroid
		for ii = 1:n_cores
			centroid_weights_split{ii} 	= calc_centroid_weights(X_split{ii}, centroid_X, parameters);
		end
	end

	% take total across difference sections of X
	centroid_weights = sum(cell2mat(centroid_weights_split), 1)';
end

%% calc_centroid_weights: calculates how many cells each centroid represents
function [centroid_weights] = calc_centroid_weights(X, centroid_X, parameters)

	% calc appropriate distance matrix
	switch parameters.metric
		case 'L1'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_L1(X, centroid_X);

		case 'L2'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_euclidean(X, centroid_X);

		case 'angle'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_angle(X, centroid_X);

		case 'corr'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_corr(X, centroid_X);

		case 'abs_corr'
			% calculate distances between points in X and points in the set of centroids ()
			D 	= comp_dist_abs_corr(X, centroid_X);

		otherwise
			error('metric not recognised')
	end

	% find minimal distances
	[min_dist, min_idx] = min(D, [], 2);
	% find where these minimal distances occur
	% min_boolean 		= bsxfun(@eq, sq_dist, min_dist);
	nn 					= size(X, 1);
	kk 					= size(centroid_X, 1);
	min_boolean 		= sparse(1:nn, min_idx, 1, nn, kk);
	% get total for each centroid
	centroid_weights 	= full(sum(min_boolean, 1));
end

%% get_cluster_labels: 
function [cluster_labels] = get_cluster_labels(X, X_split, centroids_idx, parameters)
	% fprintf('labelling each cell with closest reference node\n')

	% define useful variable
	n_cores 				= parameters.n_cores;

	% define storage variable
	cluster_labels_split 	= cell(n_cores, 1);
	centroid_X 				= X(centroids_idx, :);

	% choose parallel or not
	if parameters.pool_flag
		% find how many points in X are closest to each centroid
		parfor ii = 1:n_cores
			cluster_labels_split{ii} 	= get_cluster_labels_one(X_split{ii}, centroid_X, parameters);
		end
	else
		% find how many points in X are closest to each centroid
		for ii = 1:n_cores
			cluster_labels_split{ii} 	= get_cluster_labels_one(X_split{ii}, centroid_X, parameters);
		end
	end

	cluster_labels 	= cell2mat(cluster_labels_split);
end

%% get_cluster_labels_one: 
function [cluster_labels] = get_cluster_labels_one(X, centroid_X, parameters)
	if parameters.verbose
		show_var_size(size(X,1), 1);
	end

	% calc appropriate distance matrix
	switch parameters.metric
		case 'L1'
			% calculate distances between points in X and points in the set of centroids ()
			[~, cluster_labels] 	= comp_phi_L1(X, centroid_X);

		case 'L2'
			% calculate distances between points in X and points in the set of centroids ()
			[~, cluster_labels] 	= comp_phi_euclidean(X, centroid_X);

		case 'angle'
			% calculate distances between points in X and points in the set of centroids ()
			[~, cluster_labels] 	= comp_phi_angle(X, centroid_X);

		case 'corr'
			% calculate distances between points in X and points in the set of centroids ()
			[~, cluster_labels] 	= comp_phi_corr(X, centroid_X);

		case 'abs_corr'
			% calculate distances between points in X and points in the set of centroids ()
			[~, cluster_labels] 	= comp_phi_abs_corr(X, centroid_X);

		otherwise
			error('metric not recognised')
	end
	% % find which centroid is closest to each cluster_centroid
	% [~, cluster_labels] 	= min(D, [], 2);
end

%% comp_dist_euclidean: 
function [phi_split, idx_split] 	= comp_phi_euclidean(X, centroid_X)
	% X has size m*d, centroid_X n*d
	m 			= size(X, 1);
	n 			= size(centroid_X, 1);
	phi_split 	= zeros(m, 1);
	idx_split 	= zeros(m, 1, 'int32');

	for ii = 1:m
		% calculate squared distance to all centroids
		temp_dist_sq 	= sum((repmat(X(ii,:),n,1) - centroid_X).^2,2);

		% find minimum squared distance, return as distribution
		[phi_split(ii), idx_split(ii)] 	= min(temp_dist_sq);
	end
end

%% comp_phi_L1: 
function [phi_split, idx_split] 	= comp_phi_L1(X, centroid_X)
	% set up
	m 		= size(X, 1);
	n 		= size(centroid_X, 1);
	phi_split 	= zeros(m, 1);
	idx_split 	= zeros(m, 1, 'int32');

	for ii = 1:m
		% calculate squared distance to all centroids
		temp_dist_sq 	= sum(abs(repmat(X(ii,:), n, 1) - centroid_X),2).^2;

		% find minimum squared distance, return as distribution
		[phi_split(ii), idx_split(ii)] 	= min(temp_dist_sq);
	end
end

%% comp_phi_corr: 
function [phi_split, idx_split] 	= comp_phi_corr(X, centroid_X)
	error('corr distance not implemented yet')
	corr 	= calc_corr(X,centroid_X);
	dist 	= 1-corr; 

	% squared_distances_from_centroids 	= pdist_faster(X, centroid_X);	
	sq_dist_to_centroids 	= dist_to_centroids.^2;

	% find minimum squared distance, return as distribution
	phi_split 				= min(sq_dist_to_centroids, [], 2);
end

%% comp_phi_angle: 
function [phi_split, idx_split] 	= comp_phi_angle(X, centroid_X)
	% L2-normalization
	X 			= bsxfun(@times, X, 1./sqrt(sum(X.^2, 2)));
	centroid_X 	= bsxfun(@times, centroid_X, 1./sqrt(sum(centroid_X.^2, 2)));

	% set up
	m 			= size(X, 1);
	n 			= size(centroid_X, 1);
	phi_split 	= zeros(m, 1);
	idx_split 	= zeros(m, 1, 'int32');

	for ii = 1:m
		dot_prod 		= X(ii, :) * centroid_X';
		temp_dist_sq 	= (acos(dot_prod)/pi).^2;
		% close to zero we sometimes get complex values
		temp_dist_sq 	= real(temp_dist_sq);

		% find minimum squared distance, return as distribution
		[phi_split(ii), idx_split(ii)] 	= min(temp_dist_sq);
	end
end

%% comp_phi_abs_corr: 
function [phi_split, idx_split] 	= comp_phi_abs_corr(X, centroid_X)
	error('abs_corr distance not implemented yet')
	corr 		= calc_corr(X,centroid_X);
	dist 		= 1-abs(corr); 

	% squared_distances_from_centroids 	= pdist_faster(X, centroid_X);	
	sq_dist_to_centroids 	= dist_to_centroids.^2;

	% find minimum squared distance, return as distribution
	phi_split 				= min(sq_dist_to_centroids, [], 2);
end

%% comp_dist_euclidean: 
function dist = comp_dist_euclidean(X, centroid_X)
	% X has size m*d, centroid_X n*d
	m 		= size(X, 1);
	n 		= size(centroid_X, 1);
	dist 	= zeros(m,n);

	for ii = 1:m
		dist(ii,:) = sqrt(sum((repmat(X(ii,:),n,1) - centroid_X).^2,2)); 
	end
end

%% comp_dist_L1: 
function dist = comp_dist_L1(X, centroid_X)
	% set up
	m 		= size(X, 1);
	n 		= size(centroid_X, 1);
	dist 	= zeros(m, n);

	for ii = 1:m
		dist(ii,:) = sum(abs(repmat(X(ii,:), n, 1) - centroid_X),2); 
	end
end

%% comp_dist_corr: 
function dist = comp_dist_corr(X, centroid_X)
	corr 	= calc_corr(X,centroid_X);
	dist 	= 1-corr; 
end

%% comp_dist_angle: 
function dist = comp_dist_angle(X, centroid_X)
	% L2-normalization
	X 			= bsxfun(@times, X, 1./sqrt(sum(X.^2, 2)));
	centroid_X 	= bsxfun(@times, centroid_X, 1./sqrt(sum(centroid_X.^2, 2)));

	% dot_prod 	= X*centroid_X';
	% dist 		= acos(dot_prod)/pi;

	% set up
	m 		= size(X, 1);
	n 		= size(centroid_X, 1);
	dist 	= zeros(m, n);

	for ii = 1:m
		dot_prod 		= X(ii, :) * centroid_X';
		dist(ii, :) 	= real(acos(dot_prod)/pi);
	end
end

%% comp_dist_abs_corr: 
function dist = comp_dist_abs_corr(X, centroid_X)
	corr 		= calc_corr(X,centroid_X);
	dist 		= 1-abs(corr); 
end

%% calc_corr: 
function [corr] = calc_corr(X, centroid_X)
	% corr 	= zeros(m,n);

	% for ii = 1:m
	% 	corr(ii,:) 	= X(ii,:) * centroid_X';
	% end
	% corr 	= X*centroid_X';
	% X 			= X';
	% centroid_X 	= centroid_X';

	% zero-mean
	An 		= bsxfun(@minus, X, mean(X, 1));
	Bn 		= bsxfun(@minus, centroid_X, mean(centroid_X, 1));

	% L2-normalization
	An 		= bsxfun(@times, An, 1./sqrt(sum(An.^2, 1)));
	Bn 		= bsxfun(@times, Bn, 1./sqrt(sum(Bn.^2, 1)));

	% correlation
	corr 		= sum(An*Bn', 1);
end

%% pdist_faster: faster implementation of distance measure
function [squared_distances] = pdist_faster(X,Y)
	squared_distances = bsxfun(@plus,dot(X',X',1)',dot(Y',Y',1))-2*(X*Y');
end

%% show_var_size: display size of variables
function [] = show_var_size(m, n)
	% X has size m*d, centroid_X n*d
	size_test 	= zeros(m, n);
	size_whos 	= whos('size_test');
	clear('size_test')

	fprintf('\ndistn matrix is %d by %d, size is %.2f MB\n', m, n, size_whos.bytes / 2^20);
end
