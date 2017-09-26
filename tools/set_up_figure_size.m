%% set_up_figure_size: helper function to make figure ok for printing to pdf
function [] = set_up_figure_size(fig, units, fig_size)
	set(fig, 'units', units);
	set(fig, 'paperunits', units);
	set(fig, 'paperposition', [0, 0, fig_size]);
	set(fig, 'papersize', fig_size);
	set(fig, 'position', [0, 0, fig_size]);
	set(fig, 'paperpositionmode', 'manual');
end
