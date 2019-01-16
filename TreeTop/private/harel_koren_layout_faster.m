%% harel_koren_layout_faster: Find a nice layout for large graphs. Taken from: 
% David Harel and Yehuda Koren. A fast multi-scale method for drawing large graphs. Journal of Graph Algorithms and Applications, 2002
function [layout, energy, energy_hk] = harel_koren_layout_faster(G, options)
% function [layout] = harel_koren_layout_faster(G, options, xy)

	% admin
	check_G(G);
	if nargin == 1
		options = [];
	end 
	options 	= check_options(options);

	% get size of graph
	n_nodes 	= size(G,1);
	% Compute the all-pairs shortest path length: dV V
	shortest_path_matrix 		= dijkstra(G, (1:n_nodes)');

	% Set up a random layout L
	layout 		= rand(n_nodes, 2);
	
	% initialize neighbourhood to smallest size
	k 					= options.coarsest_graph_size;
	iteration_counter 	= 1;
	iteration_total 	= ceil((log(n_nodes) - log(options.coarsest_graph_size))/log(options.level_ratio));

	% while k ≤ |V| do
	while k <= n_nodes
		% calculate k centres across the graph
		[centres centre_idx] = calc_k_centres(shortest_path_matrix,k);
		
		% which are not centres?
		non_centres = setdiff(1:n_nodes, centres);

		% define inter-centre distance graph
		centre_shortest_paths = shortest_path_matrix(centres, centres);
		% define centre layout
		centre_layout = layout(centres, :);

		% calculate minimum distance for each centre
		nonzero_centre_shortest_paths = centre_shortest_paths;
		nonzero_centre_shortest_paths(1:k+1:k*k) = Inf;
		min_inter_centre_distances = min(nonzero_centre_shortest_paths);

		% get max of these distances (i.e. how big is largest inter-centre distance?)
		max_inter_centre_distance = max(min_inter_centre_distances);
		% radius = maxv∈centers minu∈centers{dvu} ∗ Rad
		% calculate radius as this distance times our local radius definition
		radius = max_inter_centre_distance * options.local_neighbourhood_radius;

		% LocalLayout(dcenters×centers, L(centers), radius, Iterations)
		centre_layout = calc_local_layout(centre_shortest_paths, centre_layout, radius, options.n_iterations);
		% plot_check(xy, centres, centre_layout)

		% locate each vertex at its closest centre, plus a little bit of noise (or bang on the centre for the centres themselves)
		layout(centres,:)		= centre_layout;
		layout(non_centres,:) 	= centre_layout(centre_idx(non_centres),:) + randn(numel(non_centres), 2)*options.epsilon;

		% we want the loop to be executed once for all nodes, which is why we need the following line
		if k == n_nodes
			break;
		end
		% move up to a larger neighbourhood
		k = min(k * options.level_ratio, n_nodes);
	end

	% calculate energy functions
	[energy, energy_hk] = calc_energy_values(shortest_path_matrix, layout, radius);
end

%% check_G: checks the graph 
function [] = check_G(G)
	% is it sparse?
	if ~issparse(G)
		error('G is not sparse');
	end
	% is it square
	if size(G,1) ~= size(G,2)
		error('G is not square')
	end
	% any negative values?
	neg_G = G<0;
	if any(neg_G(:))
		error('G contains negative edge weights')
	end
end

%% check_options: 
function [options] = check_options(options)
	% Define constants:

	% Rad[= 7] – determines radius of local neighborhoods
	if ~isfield(options, 'local_neighbourhood_radius')
		options.local_neighbourhood_radius = 7;
	end

	% Iterations[= 4] – determines number of iterations in local beautification
	if ~isfield(options, 'n_iterations')
		options.n_iterations = 4;
	end

	% Ratio[= 3] – ratio between number of vertices in two consecutive levels
	if ~isfield(options, 'level_ratio')
		options.level_ratio = 3;
	end

	% MinSize[= 10] – size of the coarsest graph
	if ~isfield(options, 'coarsest_graph_size')
		options.coarsest_graph_size = 10;
	end

	% define noise size around centres
	if ~isfield(options, 'epsilon')
		options.epsilon = 2e-2;
	end
end

% Goal: Find a set S ⊆ V of size k, such that max v∈V min s∈S {dsv} is minimized.
function [centres centre_idx] = calc_k_centres(shortest_path_matrix, k)
	% S �? {v} for some arbitrary v ∈ V
	S = randsample(size(shortest_path_matrix,1),1);
	% for i = 2 to k do
	for ii = 2:k	
		% 1. Find the vertex u farthest away from S
		% (i.e., such that mins∈S{dus} ≥ mins∈S{dws}, ∀w ∈ V )

		% calculate distances from S to all vertices
		all_graph_dists_to_S = shortest_path_matrix(S,:);
		if numel(S) > 1
			% find distance to closest member of S for each vertex
			[graph_dists_to_S min_idx] = min(all_graph_dists_to_S);
			% how far is furthest vertex from S?
			[max_dist_to_S furthest_vertex] = max(graph_dists_to_S);
		else
			[max_dist_to_S furthest_vertex] = max(all_graph_dists_to_S);
		end

		% 2. S �? S ∪ {u}
		% add this vertex to S
		S = unique([S furthest_vertex]);
	end

	% which is closest centre for each point?
	all_graph_dists_to_S = shortest_path_matrix(S,:);
	[graph_dists_to_S S_min_idx] = min(all_graph_dists_to_S);

	% define output values
	centres 	= S;
	centre_idx 	= S_min_idx;
end

%% initialize_delta_structure: 
function [delta_structure] = initialize_delta_structure(shortest_path_matrix, layout, in_neighbourhood_boolean)

	% define useful variables
	n_nodes = size(shortest_path_matrix,1);

	% make some matrices sparse
	in_neighbourhood_boolean 	= sparse(in_neighbourhood_boolean);
	shortest_path_matrix 		= sparse(shortest_path_matrix.*in_neighbourhood_boolean);

	% calculate x_m - x_i and y_m - y_i
	x_minus_x 			= bsxfun(@minus, repmat(layout(:,1)',n_nodes,1), repmat(layout(:,1),1,n_nodes));
	y_minus_y 			= bsxfun(@minus, repmat(layout(:,2)',n_nodes,1), repmat(layout(:,2),1,n_nodes));

	% make sparse by only considering neighbours
	x_minus_x 			= sparse(x_minus_x .* in_neighbourhood_boolean);
	y_minus_y 			= sparse(y_minus_y .* in_neighbourhood_boolean);

	% square this so we can calculate the denominator we need
	x_minus_x_squared 	= x_minus_x.^2;
	y_minus_y_squared 	= y_minus_y.^2;
	denominator_matrix	= (x_minus_x_squared + y_minus_y_squared).^0.5;
	% take inverse (but only of non-zero entries)
	[ii jj ss] 					= find(denominator_matrix);
	[mm nn] 					= size(denominator_matrix);
	inverse_denominator_matrix	= sparse(ii, jj, 1./ss, mm, nn);

	% calculate weights from distances (again restricting only to neighbours)
	[ii jj ss] 					= find(shortest_path_matrix);
	[mm nn] 					= size(shortest_path_matrix);
	weight_matrix				= sparse(ii, jj, 1./(ss.^2), mm, nn);

	% calc delta_k_v
	dE_by_dx_matrix		= weight_matrix .* (x_minus_x - shortest_path_matrix.*x_minus_x.*inverse_denominator_matrix);
	dE_by_dy_matrix		= weight_matrix .* (y_minus_y - shortest_path_matrix.*y_minus_y.*inverse_denominator_matrix);
	
	dE_by_dx 			= sum(dE_by_dx_matrix);
	dE_by_dy 			= sum(dE_by_dy_matrix);
	
	delta_k_v = full(sqrt(dE_by_dx.^2 + dE_by_dy.^2));

	% prepare output
	% permanent bits
	delta_structure.shortest_path_matrix		= shortest_path_matrix;
	delta_structure.weight_matrix				= weight_matrix;
	delta_structure.in_neighbourhood_boolean	= in_neighbourhood_boolean;

	% bits that will need updating
	delta_structure.x_minus_x					= x_minus_x;
	delta_structure.y_minus_y					= y_minus_y;
	delta_structure.inverse_denominator_matrix	= inverse_denominator_matrix;
	delta_structure.delta_k_v					= delta_k_v;
end

%% calc_x_y_adjust: 
function [adjust] = calc_x_y_adjust(delta_structure, max_idx)
	% define useful variables
	x_minus_x					= delta_structure.x_minus_x;
	y_minus_y					= delta_structure.y_minus_y;
	inverse_denominator_matrix	= delta_structure.inverse_denominator_matrix;
	shortest_path_matrix		= delta_structure.shortest_path_matrix;
	weight_matrix				= delta_structure.weight_matrix;
	in_neighbourhood_boolean	= delta_structure.in_neighbourhood_boolean;

	% calc x_m - x_i vectors
	x_minus_x_vector		= x_minus_x(:, max_idx);
	y_minus_y_vector		= y_minus_y(:, max_idx);

	% convert this into the squared denominator
	layout_x_squared 			= x_minus_x_vector.^2;
	layout_y_squared 			= y_minus_y_vector.^2;
	denominator_vector			= (layout_x_squared + layout_y_squared).^0.5;
	% take inverse (but only of non-zero entries)
	[ii jj ss] 					= find(denominator_vector);
	[mm nn] 					= size(denominator_vector);
	inverse_denominator_vector	= sparse(ii, jj, 1./ss, mm, nn);

	% define vector of shortest paths, and from this the weights
	shortest_path_vector		= shortest_path_matrix(:, max_idx);
	weight_vector 				= weight_matrix(:, max_idx);

	% calculate first derivatives of energy term
	dE_by_dx 	= sum( weight_vector .* (x_minus_x_vector - shortest_path_vector.*x_minus_x_vector	.*inverse_denominator_vector) );
	dE_by_dy 	= sum( weight_vector .* (y_minus_y_vector - shortest_path_vector.*y_minus_y_vector	.*inverse_denominator_vector) );

	% calculate second derivatives of energy term
	d2E_by_dx2 	= sum( weight_vector .* (1 - shortest_path_vector.*(y_minus_y_vector.^2)			.*(inverse_denominator_vector.^3) ) );
	d2E_by_dxdy	= sum( weight_vector .* shortest_path_vector.*x_minus_x_vector.*y_minus_y_vector	.*(inverse_denominator_vector.^3) 	);
	d2E_by_dy2 	= sum( weight_vector .* (1 - shortest_path_vector.*(x_minus_x_vector.^2)			.*(inverse_denominator_vector.^3) ) );

	% solve linear equations to get appropriate adjustment
	A = full([d2E_by_dx2 d2E_by_dxdy; d2E_by_dxdy d2E_by_dy2]);
	b = full([-dE_by_dx; -dE_by_dy]);
	adjust = (A\b)';
end

%% update_z_minus_z: 
function [out_z_minus_z] = update_z_minus_z(in_z_minus_z, max_idx, adjust)
	% [ii jj ss] 		= find(in_z_minus_z);
	% [mm nn] 		= size(in_z_minus_z);
	% update_cols		= find(ii == max_idx);
	% ss(update_cols) = ss(update_cols) - adjust;
	% update_rows		= find(jj == max_idx);
	% ss(update_rows) = ss(update_rows) + adjust;
	% out_z_minus_z	= sparse(ii, jj, ss, mm, nn);

	out_z_minus_z 				= in_z_minus_z;
	out_z_minus_z(max_idx,:)	= in_z_minus_z(max_idx,:) - adjust;
	out_z_minus_z(:,max_idx)	= in_z_minus_z(:,max_idx) + adjust;
end

%% get_chunks_from_delta_structure: 
function [x_minus_x_chunk y_minus_y_chunk inverse_denominator_matrix_chunk shortest_path_matrix_chunk weight_matrix_chunk] = get_chunks_from_delta_structure(delta_structure, neighbour_idx)
	% get chunks of these
	x_minus_x_chunk						= delta_structure.x_minus_x(:,neighbour_idx);
	y_minus_y_chunk						= delta_structure.y_minus_y(:,neighbour_idx);
	inverse_denominator_matrix_chunk	= delta_structure.inverse_denominator_matrix(:,neighbour_idx);
	shortest_path_matrix_chunk			= delta_structure.shortest_path_matrix(:,neighbour_idx);
	weight_matrix_chunk					= delta_structure.weight_matrix(:,neighbour_idx);
end

%% calc_inverse_denominator_chunk: 
function [inverse_denominator_matrix_chunk] = calc_inverse_denominator_chunk(x_minus_x_chunk, y_minus_y_chunk)
	% square this so we can calculate the denominator we need
	x_minus_x_squared 			= x_minus_x_chunk.^2;
	y_minus_y_squared 			= y_minus_y_chunk.^2;
	denominator_matrix_chunk	= (x_minus_x_squared + y_minus_y_squared).^0.5;
	% take inverse (but only of non-zero entries)
	[ii jj ss] 							= find(denominator_matrix_chunk);
	[mm nn] 							= size(denominator_matrix_chunk);
	inverse_denominator_matrix_chunk	= sparse(ii, jj, 1./ss, mm, nn);
end

%% calc_delta_k_v_chunk: 
function [delta_k_v_chunk] = calc_delta_k_v_chunk(weight_matrix_chunk, x_minus_x_chunk, y_minus_y_chunk, shortest_path_matrix_chunk, inverse_denominator_matrix_chunk)
	% calc delta_k_v
	dE_by_dx_matrix		= weight_matrix_chunk .* (x_minus_x_chunk - shortest_path_matrix_chunk.*x_minus_x_chunk.*inverse_denominator_matrix_chunk);
	dE_by_dy_matrix		= weight_matrix_chunk .* (y_minus_y_chunk - shortest_path_matrix_chunk.*y_minus_y_chunk.*inverse_denominator_matrix_chunk);
	
	dE_by_dx 			= sum(dE_by_dx_matrix);
	dE_by_dy 			= sum(dE_by_dy_matrix);
	
	delta_k_v_chunk = full(sqrt(dE_by_dx.^2 + dE_by_dy.^2));
end

%% put_chunks_back: 
function [inverse_denominator_matrix delta_k_v] = put_chunks_back(inverse_denominator_matrix, delta_k_v, inverse_denominator_matrix_chunk, delta_k_v_chunk, neighbour_idx)
	% put chunks back into big bits
	inverse_denominator_matrix(:,neighbour_idx) = inverse_denominator_matrix_chunk;
	delta_k_v(neighbour_idx) 					= delta_k_v_chunk;
end

%% prepare_delta_structure_output: 
function [delta_structure] = prepare_delta_structure_output(delta_structure, x_minus_x, y_minus_y, inverse_denominator_matrix, delta_k_v)
	% prepare output
	delta_structure.x_minus_x					= x_minus_x;
	delta_structure.y_minus_y					= y_minus_y;
	delta_structure.inverse_denominator_matrix	= inverse_denominator_matrix;
	delta_structure.delta_k_v					= delta_k_v;
end

%% udpate_delta_structure: 
function [delta_structure] = udpate_delta_structure(delta_structure, adjust, max_idx)
	% define useful variables
	in_neighbourhood_boolean	= delta_structure.in_neighbourhood_boolean;
	x_minus_x					= delta_structure.x_minus_x;
	y_minus_y					= delta_structure.y_minus_y;
	inverse_denominator_matrix	= delta_structure.inverse_denominator_matrix;
	shortest_path_matrix		= delta_structure.shortest_path_matrix;
	weight_matrix				= delta_structure.weight_matrix;
	delta_k_v					= delta_structure.delta_k_v;

	% find neighbours of max_idx
	neighbour_idx 	= in_neighbourhood_boolean(max_idx,:);

	% update x_minus_x
	x_minus_x 	= update_z_minus_z(x_minus_x, max_idx, adjust(1));
	y_minus_y 	= update_z_minus_z(y_minus_y, max_idx, adjust(2));
	
	% get chunks
	[x_minus_x_chunk y_minus_y_chunk inverse_denominator_matrix_chunk shortest_path_matrix_chunk weight_matrix_chunk] = get_chunks_from_delta_structure(delta_structure, neighbour_idx);

	% calculate update to inverse_denominator_matrix_chunk
	[inverse_denominator_matrix_chunk] = calc_inverse_denominator_chunk(x_minus_x_chunk, y_minus_y_chunk);

	% calc delta_k_v
	[delta_k_v_chunk] = calc_delta_k_v_chunk(weight_matrix_chunk, x_minus_x_chunk, y_minus_y_chunk, shortest_path_matrix_chunk, inverse_denominator_matrix_chunk);

	% put chunks back in
	[inverse_denominator_matrix delta_k_v] = put_chunks_back(inverse_denominator_matrix, delta_k_v, inverse_denominator_matrix_chunk, delta_k_v_chunk, neighbour_idx);

	% prepare output
	[delta_structure] = prepare_delta_structure_output(delta_structure, x_minus_x, y_minus_y, inverse_denominator_matrix, delta_k_v);
end

%% calc_local_layout: Find a locally nice layout L by beautifying k-neighborhoods
function [layout] = calc_local_layout(shortest_path_matrix, initial_layout, radius, n_iterations)
	% shortest_path_matrix: all-pairs shortest path length
	% initial_layout: initialized layout
	% radius: radius of neighborhoods
	
	% initialize variables
	n_nodes = size(shortest_path_matrix, 1);
	layout = initial_layout;
	in_neighbourhood_boolean = shortest_path_matrix < radius;
	% disp(['number of points not in neighbourhood = ' num2str(sum(in_neighbourhood_boolean(:)))])

	% figure;
	% hold on
	% plot(initial_layout(:,1), initial_layout(:,2), '+', 'MarkerSize', 12);

	% initialize structure we need for calculating delta_k_v and adjustments
	delta_structure = initialize_delta_structure(shortest_path_matrix, layout, in_neighbourhood_boolean);

	% for i = 1 to Iterations ∗ |V | do
	for ii = 1:( n_iterations * n_nodes )
		% 1. Choose the vertex v with the maximal ∆kv
		[max_delta max_idx] = max(delta_structure.delta_k_v);
		% 2. Compute δkv as in Kamada-Kawai
		adjust = calc_x_y_adjust(delta_structure, max_idx);
		% 3. L(v) �? L(v) + (δkv(x), δkv(y))
		layout(max_idx,:) = layout(max_idx,:) + adjust;
		% update delta values
		delta_structure = udpate_delta_structure(delta_structure, adjust, max_idx);

		% plot(layout(max_idx,1), layout(max_idx,2), '+k', 'MarkerSize', 12);
	end
	% hold off
end

%% plot_check: 
function [] = plot_check(xy, centres, centre_layout)
	xy_centres = xy(centres,:);

	labels = cellstr( num2str((1:numel(centres))') );  
	figure;
	hold on
	gscatter(xy_centres(:,1), xy_centres(:,2), centres, 'r', 'o', 24)
	text(xy_centres(:,1), xy_centres(:,2), labels)

	gscatter(centre_layout(:,1), centre_layout(:,2), centres, 'b', '+', 24)
	text(centre_layout(:,1), centre_layout(:,2), labels)

	legend('off')
	hold off
end

%% calc_energy_values: calculates how well the algorithm has done, both by standard (Kamada Kawaii) measures, and by HK's own measure
function [energy, energy_hk] = calc_energy_values(shortest_path_matrix, layout, radius)
	% calculate distances in layout
	layout_dist			= dist(layout');
	n_nodes 			= size(layout, 1);

	% calculate energies for each edge
	energy_individual 	= (layout_dist - shortest_path_matrix).^2 ./ shortest_path_matrix.^2;
	energy_individual(1:n_nodes+1:n_nodes^2) 	= 0;

	% standard energy does sum across all edges
	energy 				= sum(energy_individual(:));

	% HK energy only counts those within neighbourhoods
	is_neighbour 		= shortest_path_matrix < radius;
	energy_hk_indiv 	= energy_individual .* is_neighbour;
	energy_hk 			= sum(energy_hk_indiv(:));
end

% to do:
% - make calc_delta_k_v more efficient
% - other efficiencies?
% - check that initial layout scale is ok
% - test options



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% calc_min_graph_dist_to_S: 
% function [min_graph_dist_to_S] = calc_min_graph_dist_to_S(shortest_path_matrix, S)
% 	all_graph_dists_to_S = shortest_path_matrix(S,:);
% 	[min_graph_dist_to_S min_idx] = min(all_graph_dists_to_S);
% end

%% pdist_fastxer: faster implementation of distance measure
function [squared_distances] = pdist_faster(X,Y)
	squared_distances = bsxfun( @plus, dot(X',X',1)', dot(Y',Y',1)) - 2*(X*Y');
end

%% calc_min_dist_to_S: 
function [min_dist_to_S] = calc_min_dist_to_S(X, S)
	% calculate distances between points in X and points in the set of S
	% distances_from_S = sqrt(pdist_faster(X, S));
	distances_from_S = slmetric_pw(X', S', 'eucdist');
	% distances_from_S = pdist2(X, S).^2;

	% find minimum squared distance
	min_dist_to_S = min(distances_from_S, [], 2);
end

