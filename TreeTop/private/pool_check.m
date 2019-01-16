%% pool_check: checks that a pool is running
function [pool_flag] = pool_check(options_struct)
	if options_struct.pool_flag == false
		pool_flag 	= false;
		fprintf('running without pool\n')
		return
	end

	% is matlab version earlier than 2013b?
	version_str 	= version('-release');
	[~, idx] 		= sort({'2013b', version_str});
	old_flag 		= idx(1) == 2;

	if old_flag
		error('TreeTop requires MATLAB version 2013b or newer to run')
	else
		p 			= gcp('nocreate'); % If no pool, do not create new one.
		pool_flag 	= ~isempty(p);
	end

	if ~pool_flag
		error('TreeTop must be run with a pool (use parpool to initialize)')
	end
end
