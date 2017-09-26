%% treetop: Runs treetop for specified inputs, with specified options.
%
% Data for input into TreeTop is specified via the input_struct object. The
% required fields of input_struct are as follows:
%		data_dir		String defining path to directory where input files 
%						are stored
%		output_dir		String defining path to directory for outputs. 
%						this path has the form /parent_folder/output_name,
%						and that parent_folder already exists. The
%						output_name subdirectory is then created, and
%						output_name is used as a label for TreeTop
%						intermediate and output files.
%		used_markers	Cell array of markers to be used for this run. All 
%						must be present in the input files.
% One of the following two must be given:
%		filenames		Cell array of fcs filenames, defining input files.
%						All input files should contain the same 
%		mat_file		Inputs to TreeTop can alternatively be given via a
%						single mat file, which should contain:
% Additional possible inputs are:
%		extra_markers	Cell array of additional markers which are not used
%						in the run, but are of interest, for example for
%						validation. All must be present in the input files.
%		file_annot		Cell array of shorthand labels for fcs filenames.
%						If omitted, values in filenames are used.
%
% Options for running TreeTop are specified via the options_struct object. 
% The fields of input_struct are as follows:
%		sigma 				Bandwidth used for calculating cell density
%							values. Once TreeTop has been run, the png file
%							'[output_name] density distribution.png' can be
%							used as a diagnostic to check whether the value
%							of sigma is appropriate.
% The following options are optional
%		metric_name 		Distance metric to be used. Valid options are
%							'L1', 'L2', 'angle'; if omitted, default value
%							is 'L1'.
%		n_ref_cells  		Number of reference nodes used by TreeTop; if
%							omitted, default value is 200.
%		n_trees  			Number of trees sampled by TreeTop; if
%							omitted, default value is 1000.
%		outlier 			Quantile of density distribution below which a 
%							cell is regarded as an outlier and excluded; if
%							omitted, default value is 0.01.
%		threshold 			Quantile of density distribution above which a 
%							cell is regarded as high density, and subject 
%							to downsampling; if omitted, default value is 
%							0.5.
%		used_cofactor 		Cofactor for calculating arcsinh for 
%							used_markers. Value of 0 results in no arcsinh 
%							transformation being performed. If omitted, 
%							default value is 5.
%		extra_cofactor 		Cofactor for calculating arcsinh for 
%							extra_markers. Value of 0 results in no arcsinh 
%							transformation being performed. If omitted, 
%							default value is 5.
%		file_ext 			Specifies format of output image files. Valid
%							values are 'png', 'eps' and 'pdf'; if omitted,
%							default value is 'png'.
%
% Example code for running treetop:
%	input_struct 	= struct( ...
%		'data_dir', 			'path_to_data/wishbone clean 2', ...
%		'output_dir', 			'parent_folder_for_outputs/wishbone_clean_2_dc', ...
%		'used_markers', 		{{'DC02', 'DC03', 'DC04'}}, ...
%		'extra_markers', 		{{'CD27', 'CD4', 'CD5', 'CD127', 'CD44', 'CD69', 'CD117', 'CD62L', 'CD24', 'CD3', ...
%								'CD8', 'CD25', 'TCRb', 'BCL11b', 'CD11b', 'CD11c', 'CD161', 'CD19', 'CD38', 'CD45', ...
%								'CD90', 'Foxp3', 'GATA3', 'IA', 'Notch1', 'Notch3', 'RORg', 'Runx1', 'TCRgd', 'ki67'}}, ...
%		'filenames', 			{{'sample destiny full.fcs'}}, ...
%		'file_annot', 			{{'sample'}} ...
%		);
%	options_struct 	= struct( ...
%		'metric_name', 			'L1', ...
%		'n_ref_cells',  		200, ...
%		'n_trees',  			1000, ...
%		'outlier', 				0.01, ...
%		'threshold', 			0.2, ...
%		'sigma', 				1e-4, ...
%		'used_cofactor', 		0, ...
%		'extra_cofactor', 		0 ...
%		);
%	parpool(4)
%	treetop(input_struct, options_struct)
function treetop(input_struct, options_struct)
	fprintf('\nrunning TreeTop\n')

	% check inputs
	[input_struct, options_struct] 	= check_treetop_inputs(input_struct, options_struct);
	% check whether pool is present
	pool_check(options_struct);

	% sample ensemble of trees
	fprintf('\n1/4	Sampling ensemble of trees from data\n')
	treetop_trees(input_struct, options_struct)

	% calculate bifurcation score
	fprintf('\n2/4	Calculating branch scores for each reference node\n')
	treetop_branching_scores(input_struct, options_struct)

	% do layouts for these
	fprintf('\n3/4	Calculating layouts based on ensemble of trees\n')
	treetop_layout(input_struct, options_struct)

	% do plots
	fprintf('\n4/4	Plotting treetop outputs\n')
	treetop_plots(input_struct, options_struct)

	fprintf('\ndone.\n\n')
end
