%% save_txt_file: 
function [] = save_txt_file(save_filename, header, save_data)
	
	if isempty(header)

		dataSpec = ['%4.4f' repmat('\t%4.4f', 1, size(save_data,2) -1) '\n'];;
		fid = fopen(save_filename, 'w');
		for ii = 1:size(save_data,1)
			fprintf(fid, dataSpec, save_data(ii,:));
		end
		fclose(fid);

	else

		% if necessary, transpose header so that it has the expected dimensions
		if size(header,2) == 1
			header = header';
		elseif size(header,1) == 1

		else
			error('At least one dimension of variable header must equal one');
		end

		if size(header,2) ~= size(save_data,2)
			error(['Problem saving ' save_filename ': header and data are not compatible lengths.']);
		else
			fprintf('saving file %s\n', save_filename);
		end

		headerSpec 	= ['%s' repmat('\t%s', 1, size(header,2) -1) '\n'];
		dataSpec 	= ['%4.4f' repmat('\t%4.4f', 1, size(save_data,2) -1) '\n'];
		fid 		= fopen(save_filename, 'w');
		fprintf(fid, headerSpec, header{:});
		for ii = 1:size(save_data,1)
			fprintf(fid, dataSpec, save_data(ii,:));
		end
		fclose(fid);

	end
end

