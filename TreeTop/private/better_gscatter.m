%% better_gscatter: gscatter looks rubbish. this is a bit better
% consider adding ishold stuff
function [h_legend] = better_gscatter(x, y, g, options)
	% sort out holding
	hold_status 	= ishold;
	if ~hold_status
		hold on
	end

	% get palette
	ctype 		= 'qual';
	palette 	= 'Set1';
	point_size 	= 10;
	legend_flag = false;
	location 	= 'EastOutside';

	if nargin > 3 
		if isfield(options, 'palette')
			palette 	= options.palette;
		end
		if isfield(options, 'size')
			point_size 	= options.size;
		end
		if isfield(options, 'legend')
			legend_flag = options.legend;
		end
		if isfield(options, 'leg_pos')
			location 	= options.leg_pos;
		end
	end

	% get names
	[g_vals, ~, g_idx] 	= unique(g);
	if ~iscell(g_vals)
		g_vals 		= arrayfun(@num2str, g_vals, 'unif', false);
	end

	% get g
	n_cols 			= length(g_vals);
	h 				= zeros(n_cols, 1);

	if n_cols > 11
		max_cols 		= 11;
		pal_mat_short 	= cbrewer('div', 'RdYlBu', max_cols);
		x_vals 			= (0:max_cols - 1) / (max_cols - 1);
		x_queries 		= (0:n_cols - 1) / (n_cols - 1);
		pal_mat 		= arrayfun(@(jj) interp1(x_vals, pal_mat_short(:,jj), x_queries)', 1:size(pal_mat_short,2), 'uniformoutput', false);
		pal_mat 		= [pal_mat{:}];

	elseif n_cols > 9 & n_cols <= 11
		pal_mat 		= cbrewer('div', 'RdYlBu', n_cols);
		
	elseif n_cols < 3
		pal_mat 		= cbrewer(ctype, palette, 3);
		
	else
		pal_mat 		= cbrewer(ctype, palette, n_cols);
		
	end

	for ii = 1:n_cols
		this_idx 	= g_idx == ii;
		h(ii) 		= plot(x(this_idx), y(this_idx), '.', 'color', pal_mat(ii, :), 'markersize', point_size*2);
	end

	if legend_flag
		h_legend 		= legend(h, g_vals{:}, 'Location', location);
	else
		h_legend 		= [];
	end

	if ~hold_status
		hold off
	end
end
