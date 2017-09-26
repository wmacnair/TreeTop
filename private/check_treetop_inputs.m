%% CHECK_TREETOP_INPUTS: Checks that inputs to TreeTop are ok. Adds default values where none given.
% [input_struct, options_struct] = CHECK_TREETOP_INPUTS(input_struct, options_struct) checks both input_struct
% and options_struct, and amends them with default values where necessary.
% [input_struct, ~] = CHECK_TREETOP_INPUTS(input_struct, []) checks only input_struct, returning
% an empty object for options_struct.
% [~, options_struct] = CHECK_TREETOP_INPUTS([], options_struct) checks only options_struct,
% returning an empty object for input_struct.
function [input_struct, options_struct] = check_treetop_inputs(input_struct, options_struct)
	if exist('input_struct', 'var')
		input_struct 	= check_input_struct(input_struct);
	else
		fprintf('no input_struct given as input; returning empty value\n')
		input_struct 	= [];
	end

	if exist('options_struct', 'var')
		options_struct 	= check_options_struct(options_struct);
	else
		fprintf('no options_struct given as input; returning empty value\n')
		options_struct 	= [];
	end
end

%% check_input_struct: 
function [input_struct] = check_input_struct(input_struct)
	if isempty(input_struct)
		error('input_struct is empty; check that extra_markers is specified as {{}}, not {}')
	end
	% check all fieldnames
	fields_list 	= fieldnames(input_struct);
	input_required 	= {'data_dir', 'output_dir', 'used_markers'};
	missing_fields 	= setdiff(input_required, fields_list);
	if ~isempty(missing_fields)
		error(sprintf('the following fields are missing from input_struct: %s', strjoin(missing_fields, ', ')))
	end
	if ~isfield(input_struct, 'filenames') & ~isfield(input_struct, 'mat_file')
		error('input_struct must have one of mat_file or filenames as field names')
	end

	% 	% save_stem can't have branch in it
	% 	if ~isempty(regexp(input_struct.save_stem, '^branch'))
	% 		error('save_stem can''t start with "branch"')
	% 	end

	% check data_dir
	if ~exist(input_struct.data_dir, 'dir')
		error('data_dir does not exist\n')
	end
	% if fcs files specified, check that they exist
	
	% if mat_file specified, check that it exists
	if isfield(input_struct, 'mat_file') && ~exist(fullfile(input_struct.data_dir, input_struct.mat_file), 'file')
		error('specified mat_file does not exist')
	end

	% add save_stem from output directory
	[parent_dir, name, ext] 	= fileparts(input_struct.output_dir);
	% deal with full stops
	if ~isempty(ext)
		name 	= [name, ext];
	end
	input_struct.save_stem 	= name;

	% check parent_dir
	if ~exist(parent_dir, 'dir')
		error('parent directory of output_dir does not exist\n')
	end

	% make output_dir
	if ~exist(input_struct.output_dir, 'dir')
		mkdir(input_struct.output_dir)
	end

	% if extra_markers missing, make it empty
	if ~isfield(input_struct, 'extra_markers')
		input_struct.extra_markers 	= {};
	elseif any(cellfun(@isempty, input_struct.extra_markers))
		error('some extra_markers are empty (to have no extra markers, use {}, or leave extra_markers undefined')
	end

	% if file_annot missing, make it same as filenames
	if isfield(input_struct, 'filenames') & ~isfield(input_struct, 'file_annot')
		input_struct.file_annot 	= {input_struct.filenames};
	end

	% if cofactor missing, make it 5
	if ~isfield(input_struct, 'used_cofactor')
		input_struct.used_cofactor 	= 5;
	end
	if ~isfield(input_struct, 'extra_cofactor')
		input_struct.extra_cofactor = 5;
	end
end

%% check_options_struct: 
function [options_struct] 	= check_options_struct(options_struct);
	% same for options:
	fields_list 		= fieldnames(options_struct);
	options_required 	= {'sigma'};
	missing_fields 		= setdiff(options_required, fields_list);
	if ~isempty(missing_fields)
		error(sprintf('the following fields are missing from options_struct: %s', strjoin(missing_fields, ', ')))
	end

	% set defaults if missing
	if ~isfield(options_struct, 'metric_name')
		options_struct.metric_name 	= 'L1';
	end
	if ~isfield(options_struct, 'outlier')
		options_struct.outlier		= 0.01;
	end
	if ~isfield(options_struct, 'threshold')
		options_struct.threshold 	= 0.5;
	end
	if ~isfield(options_struct, 'sample_size')
		options_struct.sample_size 	= 1e5;
	end
	if ~isfield(options_struct, 'n_ref_cells')
		options_struct.n_ref_cells 	= 200;
	end
	if ~isfield(options_struct, 'n_dens_ref')
		options_struct.n_dens_ref	= 5000;
	end
	if ~isfield(options_struct, 'n_trees')
		options_struct.n_trees 		= 1000;
	end
	if ~isfield(options_struct, 'file_ext')
		options_struct.file_ext 	= 'png';
	end
	if ~isfield(options_struct, 'seed')
		options_struct.seed 		= 73;
	end
	if ~isfield(options_struct, 'p_cutoff')
		options_struct.p_cutoff 	= 0.95;
	end
	if ~isfield(options_struct, 'three_flag')
		options_struct.three_flag 	= false;
	end
	if ~isfield(options_struct, 'pool_flag')
		options_struct.pool_flag 	= true;
	end
end
