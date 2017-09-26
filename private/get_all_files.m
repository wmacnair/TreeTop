%% get_all_files: loads all files
function all_struct = get_all_files(input_struct)

	% either get from mat file
	if isfield(input_struct, 'mat_file')
		all_struct 		= get_from_mat_file(input_struct);
		return
	end
	% or from fcs

	% unpack 
	filenames 			= input_struct.filenames;
	file_annot 			= input_struct.file_annot;
	used_markers 		= input_struct.used_markers;
	extra_markers 		= input_struct.extra_markers;

	% define storage variables
	n_files 			= length(input_struct.filenames);
	used_cell 			= cell(n_files, 1);
	extra_cell 			= cell(n_files, 1);
	celltype_cell 		= cell(n_files, 1);

	% celltype file for index
	celltype_file 		= fullfile(input_struct.data_dir, filenames{1});
	[~, celltype_hdr] 	= fca_readfcs_3_1(celltype_file);
	marker_names 		= {celltype_hdr.par.name2}';
	empty_idx			= cellfun(@isempty, marker_names);
	non_empty_markers	= marker_names(~empty_idx);
	
	% if name2 didn't work, then try name
	if numel(intersect(used_markers, non_empty_markers)) ~= numel(used_markers)
		marker_names    	= {celltype_hdr.par.name}';
		empty_idx			= cellfun(@isempty, marker_names);
		non_empty_markers	= marker_names(~empty_idx);
	end

	% check we can find all requested markers in the first file
	if numel(intersect(used_markers, non_empty_markers)) ~= numel(used_markers)
		error('file marker names don''t match those requested in input_struct')
	end
	
	% match used and extra markers to fcs files, check no duplicate marker names
	used_idx 			= cellfun(@(test_str) find(strcmp(test_str, marker_names)), used_markers, 'unif', false);
	if ~all(cellfun(@length, used_idx) == 1)
		error('some names in used_markers match fcs file more than once')
	end
	extra_idx 			= cellfun(@(test_str) find(strcmp(test_str, marker_names)), extra_markers, 'unif', false);
	if ~all(cellfun(@length, extra_idx) == 1)
		error('some names in extra_markers match fcs file more than once')
	end

	% check we can find all the markers we requested
	used_idx 			= cell2mat(used_idx);
	extra_idx 			= cell2mat(extra_idx);
	if numel(used_idx) ~= numel(used_markers) | numel(extra_idx) ~= numel(extra_markers)
		error('markers didn''t match those in file...')
	end

	% loop through all files
	fprintf('opening %d files:\n', n_files)
	for ii = 1:n_files
		fprintf('.')

		this_file 			= fullfile(input_struct.data_dir, filenames{ii});
		[fcsdat, fcshdr] 	= fca_readfcs_3_1(this_file);

		% restrict to markers of interest
		used_data 			= fcsdat(:, used_idx);
		extra_data 			= fcsdat(:, extra_idx);

		% do arcsinh
		used_data 			= flow_arcsinh(used_data, input_struct.used_cofactor);
		extra_data 			= flow_arcsinh(extra_data, input_struct.extra_cofactor);
		n_cells 			= size(used_data, 1);

		% store
		used_cell{ii} 		= used_data;
		extra_cell{ii} 		= extra_data;
		celltype_cell{ii} 	= repmat(file_annot(ii), n_cells, 1);
	end
	fprintf('\n')

	% put all together
	fprintf('combining into one matrix\n')
	used_data 			= cell2mat(used_cell);
	extra_data 			= cell2mat(extra_cell);
	all_celltype 		= vertcat(celltype_cell{:});

	% assemble into output struct
	all_struct 	= struct( ...
		'used_data', 		{used_data}, ...
		'extra_data', 		{extra_data}, ...
		'all_labels', 		{all_celltype}, ...
		'used_markers', 	{used_markers}, ...
		'extra_markers', 	{extra_markers} ...
		);
end

%% get_from_mat_file: 
function all_struct = get_from_mat_file(input_struct)
	% unpack
	used_markers	= input_struct.used_markers;
	extra_markers	= input_struct.extra_markers;
	
	% load file
	load(fullfile(input_struct.data_dir, input_struct.mat_file), 'all_struct')

	% double check a couple of things
	fields_list 	= fieldnames(all_struct);
	input_required 	= {'all_data', 'all_markers', 'all_labels'};
	missing_fields 	= setdiff(input_required, fields_list);
	if ~isempty(missing_fields)
		error(sprintf('the following fields are missing from the variable all_struct in file %s:\n%s', mat_file, strjoin(missing_fields, ', ')))
	end

	% then subdivide appropriately
	all_markers 	= all_struct.all_markers;
	used_idx 		= ismember(all_markers, input_struct.used_markers);
	n_missing		= length(used_markers) - sum(used_idx);
	if n_missing > 0
		error('%d markers not found in data', n_missing)
	end
	all_struct 					= all_struct;
	all_struct.used_data 		= all_struct.all_data(:, used_idx);
	all_struct.used_markers  	= input_struct.used_markers;
	all_struct.all_labels  		= all_struct.all_labels;

	if isempty(extra_markers)
		all_struct.extra_data  		= [];
		all_struct.extra_markers  	= {};
	else
		extra_idx					= ismember(all_markers, input_struct.extra_markers);
		if sum(extra_idx) < length(extra_markers)
			error('missing extra marker')
		end
		if sum(extra_idx) > length(extra_markers)
			error('extra extra marker')
		end
		all_struct.extra_data  		= all_struct.all_data(:, extra_idx);
		all_struct.extra_markers  	= input_struct.extra_markers;

		n_missing		= length(used_markers) - sum(used_idx);
		if n_missing > 0
			error('%d markers not found in data', n_missing)
		end
	end
end

%% flow_arcsinh: 
% cofactor can be one signal number, for all channels
%          can be a vector, one element for one channel
% if cofactor == 0, then it means linear, no transform. 
function [data] = flow_arcsinh(data, cofactor)

	if length(cofactor)~=1 && length(cofactor) ~= size(data, 2)
		error('wrong input of cofactors; data not transformed')
	end

	if length(cofactor)==1
		if cofactor==0
			return
		else
			data = log( data(:,:)/cofactor + sqrt((data(:,:)/cofactor).^2+1) );
			return
		end
	end

	if length(cofactor) == size(data,1)
		for i=1:size(data,1)
			if cofactor(i)==0
				continue;
			else
				data(i,:) = log( data(i,:)/cofactor(i) + sqrt((data(i,:)/cofactor(i)).^2+1) );
			end        
		end
	end
end
