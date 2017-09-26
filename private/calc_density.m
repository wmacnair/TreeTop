%% calc_density: calculates density values
function [density_vector] = calc_density(sample_struct, options_struct)
	% do we want to save some of these outputs?
	kmedoids_flag 		= isfield(options_struct, 'kmedoids_flag') && options_struct.kmedoids_flag;
	% define how big a reference to do
	% and sigma values

	% unpack
	used_data 			= sample_struct.used_data;
	sigma 				= options_struct.sigma;
	metric_name 		= options_struct.metric_name;
	n_dens_ref 			= options_struct.n_dens_ref;
	n_samples			= size(used_data, 1);
	n_dens_ref 			= min(n_dens_ref, n_samples);

	% define maximum matrix size
	max_mat_entries 	= 1e6;

	% split data into chunks
	chunk_cell 			= calculate_chunks(used_data, n_dens_ref, max_mat_entries);
	chunk_idx 			= cell2mat(chunk_cell);

	% sample reference cells for calculating density
	sample_idx 			= randsample(n_samples, n_dens_ref);
	dens_ref_mat 		= used_data(sample_idx, :);

	% loop through chunks
	n_chunks 			= numel(chunk_cell);
	density_cell 		= cell(n_chunks, 1);

	% CHECK THIS
	kk 					= 100;
	if kmedoids_flag
		knn_chunk 		= cell(n_chunks, 1);
	end	

	fprintf('calculating density for all points in sample\n')
	if options_struct.pool_flag
		% define function for each chunk
		parfor ii = 1:n_chunks
			% restrict samples to just this chunk
			this_idx 			= find(chunk_idx == ii);
			this_mat 			= used_data(this_idx, :);

			% calculate distance matrix
			dist_mat 			= all_distance_fn(this_mat, dens_ref_mat, options_struct.metric_name);
			% zero_idx			= dist_mat(:) == 0;
			% dist_mat(zero_idx)	= Inf;
			
			% check whether they overlap with dens sample
			overlap_idx 		= ismember(this_idx, sample_idx);
			% if sum(overlap_idx) ~= sum(zero_idx)
			% 	error('overlaps don''t match')
			% end

			% do gaussian kernel of these for each sigma
			gauss_dist 			= sum(exp( -(dist_mat/sigma).^2 /2 ), 2);
			% remove value of one where there's an overlap
			gauss_dist 			= gauss_dist - overlap_idx;

			% store
			density_cell{ii} 	= gauss_dist;

			if kmedoids_flag
				% also do knn density
				knn_dens 		= arrayfun(@(ii) kk_th_val(dist_mat(ii, :), kk), 1:size(dist_mat, 1));
				knn_chunk{ii} 	= knn_dens;
			end
		end
	else
		% define function for each chunk
		for ii = 1:n_chunks
			% restrict samples to just this chunk
			this_idx 			= find(chunk_idx == ii);
			this_mat 			= used_data(this_idx, :);

			% calculate distance matrix
			dist_mat 			= all_distance_fn(this_mat, dens_ref_mat, options_struct.metric_name);
			% zero_idx			= dist_mat(:) == 0;
			% dist_mat(zero_idx)	= Inf;
			
			% check whether they overlap with dens sample
			overlap_idx 		= ismember(this_idx, sample_idx);
			% if sum(overlap_idx) ~= sum(zero_idx)
			% 	error('overlaps don''t match')
			% end

			% do gaussian kernel of these for each sigma
			gauss_dist 			= sum(exp( -(dist_mat/sigma).^2 /2 ), 2);
			% remove value of one where there's an overlap
			gauss_dist 			= gauss_dist - overlap_idx;

			% store
			density_cell{ii} 	= gauss_dist;

			if kmedoids_flag
				% also do knn density
				knn_dens 		= arrayfun(@(ii) kk_th_val(dist_mat(ii, :), kk), 1:size(dist_mat, 1));
				knn_chunk{ii} 	= knn_dens;
			end
		end
	end

	% put into vector
	density_vector 		= cell2mat(density_cell);

	% save knn density
	if kmedoids_flag
		% put all knns together
		knn_dens 		= cell2mat(knn_chunk);

		% save outputs
		save(fullfile(options_struct.outlier_dir, 'dens_values.mat'), 'density_vector', 'knn_dens')
	end

	% check sizes ok
	if size(density_vector, 1) ~= n_samples
		error('density_vector doesn''t have same number of entries as used_data')
	end
	% check no NaNs
	if any(isnan(density_vector))
		error('NaN in density_vector')
	end
end

%% calculate_chunks: 
function chunk_cell = calculate_chunks(used_data, n_dens_ref, max_mat_entries)
	% how big should chunks be?
	n_samples 			= size(used_data, 1);
	n_dens_ref 			= min(n_dens_ref, n_samples);
	chunk_size			= floor(max_mat_entries / n_dens_ref);

	% want at most M entries in matrix
	% n_dens_ref * chunk_size <= M
	% so divide into chunks of size at most 
	n_chunks 			= ceil(n_samples / chunk_size);
	remainder 			= n_chunks * chunk_size - n_samples;
	size_vector 		= repmat(chunk_size, n_chunks, 1);
	size_vector(end) 	= chunk_size - remainder;
	% double check that size_vector has right number
	if sum(size_vector) ~= n_samples
		error('chunks wrong')
	end
	% make chunk indices
	chunk_cell 			= arrayfun(@(ii) repmat(ii, size_vector(ii), 1), (1:n_chunks)', 'unif', false);
	chunk_idx 			= cell2mat(chunk_cell);
	if numel(chunk_idx) ~= n_samples
		error('chunks wrong')
	end
end

%% kk_th_val: 
function [val] = kk_th_val(row, kk)
	sorted_row 	= sort(row);
	val 		= sorted_row(kk);
end