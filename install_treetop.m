%% install_treetop: function to add treetop to the path, and check that several necessary compiled functions are functional
function [] = install_treetop()
	% set up path
	treetop_path()

	% check compiled functions
	check_compiled_functions()
end

function treetop_path()
	% save changes to path and startup.m
	prompt 	= 'TreeTop needs to be added to MATLAB''s path to run.\nWould you like to save this change to the path?\n(Saving the change means the path variable is slightly larger.\nNot saving means that you will need to run this function again, the next time you want to use TreeTop.)\nY/N [default is N]:';
	str 	= input(prompt, 's');
	if isempty(str)
		str = 'N';
	end
	str 	= upper(str);
	if ~ismember(str, {'Y', 'N'})
		error('input must be one of: Y, y, N, n')
	end
	SAVE_PERMANENTLY 	= strcmp(str, 'Y');
	
	% add paths
	[install_dir, ~, ~] = fileparts(mfilename('fullpath'));
	addpath(genpath(install_dir));
	if (SAVE_PERMANENTLY)
		savepath;
	end
end

function check_compiled_functions()
	% check this function
	dijk_mex 	= check_mex_function('dijkstra.cpp', './TreeTop/private');
	mi_mex 		= check_mi_toolbox('MIToolboxMex.c', './MIToolbox/matlab');

	% check each function works
	dijk_works 	= false;
	if dijk_mex
		dijk_works 	= check_dijkstra_works();
	end
	mi_works 	= false;
	if mi_mex
		mi_works 	= check_mi_works();
	end

	% report results
	if dijk_works & mi_works
		fprintf('All compiled files available and working; TreeTop should run without any problems\n');
	elseif dijk_works & ~mi_works
		fprintf('The relevant compiled file for the mutual information toolbox was either not available or not working.\n')
		fprintf('treetop should run ok, but outputs from the optional treetop_pre_run script will be incomplete.\n')
		fprintf('To fix this, make sure you have a valid compiler installed and\n')
		fprintf('compile the function ./MIToolbox/matlab/MIToolboxMex.\n');
	elseif ~dijk_works
		fprintf('The relevant compiled file for the function dijkstra was either not available or not working.\n')
		fprintf('treetop needs this to run; to fix this, make sure you have a valid compiler installed and\n')
		fprintf('compile the function ./TreeTop/private/dijkstra.\n');
	end
end

function mex_exists = check_mex_function(fn_name, fn_dir)
	start_dir 	= cd;
	cd(fn_dir)
	fprintf('testing whether function %s is compiled\n', fn_name)

	% Check input filename
	assert(ischar(fn_name),'source_file must be a string')

	% Check extension is specified
	assert(~isempty(strfind(fn_name,'.')),'source_file: no file extension specified')

	% Locate source file
	[pathstr, name, ext] 	= fileparts(which(fn_name));
	% Create filename, mex filename
	filename 	= [pathstr filesep name ext];
	mexfilename = [pathstr filesep name '.' mexext];

	% check whether source exists
    mex_exists  = false;
	if strcmp(pathstr,'')
		% source file not found
		fprintf([source_file ': not found'])

	% now check whether mexfile exists
	elseif exist(mexfilename, 'file') ~= 3
		 % if source file does not exist, try to compile
		disp(['Compiling "' name ext '".'])

		% compile, with options if appropriate
		try
			mex(fn_name)
			fprintf('Function %s successfully compiled\n', fn_name);
			mex_exists 	= true;

		catch lasterr
			fprintf('Could not compile function %s. \n', fn_name);
        end		
    else
        mex_exists  = true;
    end

	% switch back to original directory
	cd(start_dir)
end

function mex_exists = check_mi_toolbox(fn_name, fn_dir)
	start_dir 	= cd;
	cd(fn_dir)
	fprintf('testing whether function %s is compiled\n', fn_name)

	% Check input filename
	assert(ischar(fn_name),'source_file must be a string')

	% Check extension is specified
	assert(~isempty(strfind(fn_name,'.')),'source_file: no file extension specified')

	% Locate source file
	[pathstr, name, ext] 	= fileparts(which(fn_name));
	% Create filename, mex filename
	filename 	= [pathstr filesep name ext];
	mexfilename = [pathstr filesep name '.' mexext];

	% check whether source exists
    mex_exists  = false;
	if strcmp(pathstr,'')
		% source file not found
		fprintf([source_file ': not found'])

	% now check whether mexfile exists
	elseif exist(mexfilename, 'file') ~= 3
		 % if source file does not exist, try to compile
		disp(['Compiling "' name ext '".'])

		% compile, with options if appropriate
		try
			CompileMIToolbox
			fprintf('Function %s successfully compiled\n', fn_name);
			mex_exists 	= true;

		catch lasterr
			fprintf('Could not compile function %s. \n', fn_name);
        end		
    else
        mex_exists  = true;
    end

	% switch back to original directory
	cd(start_dir)
end

%% check_dijkstra_works: 
function [dijk_works] = check_dijkstra_works()
	dijk_works 	= false;

	% define things to test with
	n_nodes 	= 10;
	G 			= zeros(n_nodes);
	for ii = 1:n_nodes-1
		G(ii, ii+1) 	= 1;
		G(ii+1, ii) 	= 1;
	end
	G 			= sparse(G);
	D_true 		= zeros(n_nodes);
	for ii = 1:n_nodes
		D_true(ii, :) 	= abs((1:n_nodes) - ii);
	end

	% check dijkstra works on simple test case
	current_dir = cd;
	cd('./TreeTop/private')
	try
		D 			= dijkstra(G, (1:n_nodes)');
		dijk_works 	= all(all(abs(D - D_true) < 1e-14));
	catch lasterr
		% fprintf('dijkstra did not work\n')
	end
	cd(current_dir)
end

%% check_mi_works: 
function [mi_works] = check_mi_works()
	mi_works 	= false;

	% define things to test with
	n_entries 	= 10;
	X 			= repmat(1:n_entries, 1, n_entries);
	mi_true 	= n_entries;

	% check MI works on simple test case
	current_dir = cd;
	cd('./TreeTop/private')
	try
		mi_val 		= MIToolboxMex(7, X', X');
		mi_works 	= abs(mi_true - 2^mi_val) < 1e-14;
	catch lasterr
		fprintf('MI toolbox did not work\n')
	end
	cd(current_dir)
end

% Max's comments:
% Firstly: I would add a disclaimer at the start that matlab needs to have a C compiler installed (which has to be done manually for mac)
% Secondly: regarding the output_dir variable which holds the string for the output path. You state that this path "has the form /parent_folder/output_name," with parent_folder existing but output_name not and only being generated by the pre run function. For me the whole path needs to already exist otherwise it throws an error.
% Regarding the treetop_pre_run function: here I first had to go to: path/TreeTop-master/MIToolbox/matlab and manually compile the MIToolbox.m function otherwise it would throw an error. (It does not compile automatically)
% Then when running treetop(input_struct, options_struct): I had to rename the (TreeTop-master/TreeTop/private folder because I could not add it to the path when it was named private^^
% Also the compiled version of the dijkstra.m inside TreeTop/private does not work.
% I had to download the whole Toolbox and then also manually compile the dijkstra function. Then I had to make sure that treetop uses this version and not the one inside /TreeTop/private
% But thats already all. Everything else runs smoothly

