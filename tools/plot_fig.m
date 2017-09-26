%% plot_fig: 
function plot_fig(fig, name_stem, file_ext, fig_size)
	plot_file 		= sprintf('%s.%s', name_stem, file_ext);

	% set up figure
	units 			= 'inches';
	set_up_figure_size(fig, units, fig_size)

	% do plot
	switch file_ext
		case 'png'
			r 			= 300; % pixels per inch
			print(fig, '-dpng', sprintf('-r%d', r), plot_file);
		case 'eps'
			print(fig, '-depsc', plot_file);
		case 'pdf'
			print(fig, '-dpdf', plot_file);
		otherwise
			error('invalid plot extension')
	end
	close(fig)
end

%% set_up_figure_size: helper function to make figure ok for printing to pdf
function [] = set_up_figure_size(fig, units, fig_size)
	set(fig, 'units', units);
	set(fig, 'paperunits', units);
	set(fig, 'paperposition', [0, 0, fig_size]);
	set(fig, 'papersize', fig_size);
	set(fig, 'position', [0, 0, fig_size]);
	set(fig, 'paperpositionmode', 'manual');
end
