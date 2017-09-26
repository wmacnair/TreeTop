%% pool_check: checks that a pool is running
function [pool_flag] = pool_check(options_struct)
	if options_struct.pool_flag == false
		pool_flag 	= false;
		fprintf('running without pool\n')
		return
	end

	% is matlab version earlier than 2013b?
	warning('should we limit which versions it can run in?')
	version_str 	= version('-release');
	[~, idx] 		= sort({'2013b', version_str});
	old_flag 		= idx(1) == 2;

	if old_flag
		num_workers = matlabpool('size');
		pool_flag 	= num_workers > 1;
	else
		p 			= gcp('nocreate'); % If no pool, do not create new one.
		pool_flag 	= ~isempty(p);
	end

	if ~pool_flag
		error('treetop must be run with a pool (use parpool to initialize)')
	end
end
