%% treetop_example_runs: Set of runs which reproduce the results in the paper.
% These are based on the entries in zip file treetop_data. 
% data_dir is the parent path of where treetop_data was unzipped.
% output_dir is the parent path of where you would like outputs to be stored.
% run_switch is the index of which run you would like to do, out of the following:
%   1 	T cell thymic maturation data, preprocessed with diffusion maps
%   2 	Healthy human bone marrow, mass cytometry data
%   3 	Hierarchically branching synthetic data (generated by Stefan Ganscha)
%   4 	Linear synthetic data
%   5 	Gaussian synthetic data
%   6 	Circular synthetic data
%   7 	Triangular synthetic data
%   8 	B cell maturation
%   9 	Healthy human bone marrow, scRNAseq data (Paul et al., 2015)
%   10 	Healthy human bone marrow, scRNAseq data (Velten et al., 2017)
%   11 	Swiss roll synthetic data
function treetop_example_runs(data_dir, output_dir, run_switch)
	% set up pool
	current_pool 	= gcp('nocreate');
	if numel(current_pool) == 0
		error('Please call parpool before running treetop_example_runs')
	end

	switch run_switch			
		case 1
			% T cell thymic maturation data, preprocessed with diffusion maps
			% Setty et al. 2016. “Wishbone Identifies Bifurcating Developmental Trajectories from Single-Cell Data.” Nature Biotechnology, May. doi:10.1038/nbt.3569.
			fprintf('running treetop on T cell thymic maturation data, preprocessed with diffusion maps\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'thymus'), ...
				'output_dir', 			fullfile(output_dir, 'thymus'), ...
				'used_markers', 		{{'DC02', 'DC03', 'DC04'}}, ...
				'extra_markers', 		{{'CD27', 'CD4', 'CD5', 'CD127', 'CD44', 'CD69', 'CD117', 'CD62L', 'CD24', 'CD3', ...
										'CD8', 'CD25', 'TCRb', 'BCL11b', 'CD11b', 'CD11c', 'CD161', 'CD19', 'CD38', 'CD45', ...
										'CD90', 'Foxp3', 'GATA3', 'IA', 'Notch1', 'Notch3', 'RORg', 'Runx1', 'TCRgd', 'ki67'}}, ...
				'filenames', 			{{'wishbone thymus 2.fcs'}}, ...
				'file_annot', 			{{'thymus'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.2, ...
				'sigma', 				1e-4, ...
				'layout_seed', 			2 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 2
			% Healthy human bone marrow, mass cytometry data
			% Amir et al. 2013. “viSNE Enables Visualization of High Dimensional Single-Cell Data and Reveals Phenotypic Heterogeneity of Leukemia.” Nature Biotechnology 31 (6).
			fprintf('running treetop on healthy human bone marrow, mass cytometry data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'bone marrow cytof'), ...
				'output_dir', 			fullfile(output_dir, 'bone marrow cytof'), ...
				'used_markers', 		{{'CD11c(Tb159)Di', 'CD14(Gd160)Di', 'CD33(Yb173)Di', 'CD3(Er170)Di', 'CD45(Sm154)Di', 'CD15(Dy164)Di', ...
										'CD24(Dy161)Di', 'CD19(Nd142)Di', 'CD22(Nd143)Di', 'CD20(Sm147)Di', 'CD117(Yb171)Di', 'IgM-s(Lu175)Di', ...
										'IgM-i(Eu153)Di', 'HLADR(Yb174)Di', 'CD79b(Nd146)Di', 'CD38(Er168)Di', 'CD235-62-66b(In113)Di', ...
										'CD72(Eu151)Di', 'CD7(Yb176)Di', 'CD47(Nd145)Di'}}, ...
				'extra_markers', 		{{'CD49d(Yb172)Di', 'Pax5(Ho165)Di', 'CD127(Dy162)Di', 'TdT(Dy163)Di', 'CD34(Nd148)Di', 'CD10(Gd156)Di', ...
										'CD179b(Gd158)Di', 'CD179a(Sm149)Di'}}, ...
				'filenames', 			{{'bone_marrow_T_cells.fcs', 'bone_marrow_Myeloid.fcs', ...
										'bone_marrow_ungated_CD7_low.fcs', 'bone_marrow_CD24_hi.fcs', ...
										'bone_marrow_B_cells.fcs', 'bone_marrow_NK_cells.fcs', ...
										'bone_marrow_HSCs.fcs', 'bone_marrow_ungated_CD20hi_-_CD3hi.fcs'}}, ...
				'file_annot', 			{{'T cells', 'Myeloid', 'ungated CD7 low', 'Granulocytes', 'B cells', 'NK cells', 'HSCs', 'ungated CD20hi - CD3hi'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				2 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 3
			% Hierarchically branching synthetic data (generated by Stefan Ganscha)
			% Ocone et al. 2015. “Reconstructing Gene Regulatory Dynamics from High-Dimensional Single-Cell Snapshot Data.” Bioinformatics  31 (12): i89–96.
			fprintf('running treetop on hierarchically branching synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic hierarchical'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic hierarchical'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('D%02d', ii), 1:12, 'unif', false)}, ...
				'extra_markers', 		{{'time'}}, ...
				'filenames', 			{{'synthetic hierarchically branching.fcs'}}, ...
				'file_annot', 			{{'branching'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				1 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 4
			% Linear synthetic data
			fprintf('running treetop on linear synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic linear'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic linear'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('D%02d', ii), 1:10, 'unif', false)}, ...
				'used_cofactor', 		0, ...
				'filenames', 			{{'synthetic linear.fcs'}}, ...
				'file_annot', 			{{'linear'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				0.05 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)

		case 5
			% Gaussian synthetic data
			fprintf('running treetop on Gaussian synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic gaussian'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic gaussian'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('D%02d', ii), 1:10, 'unif', false)}, ...
				'used_cofactor', 		0, ...
				'filenames', 			{{'synthetic gaussian.fcs'}}, ...
				'file_annot', 			{{'gaussian'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				1 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)

		case 6
			% Circular synthetic data
			fprintf('running treetop on circular synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic circle'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic circle'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('readout%d', ii), 1:10, 'unif', false)}, ...
				'extra_markers', 		{{'angle'}}, ...
				'used_cofactor', 		0, ...
				'extra_cofactor', 		0, ...
				'filenames', 			{{'synthetic circle.fcs'}}, ...
				'file_annot', 			{{'circle'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				0.05 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)

		case 7
			% Triangular synthetic data
			fprintf('running treetop on triangular synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic triangle'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic triangle'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('D%02d', ii), 1:10, 'unif', false)}, ...
				'used_cofactor', 		0, ...
				'filenames', 			{{'synthetic triangle.fcs'}}, ...
				'file_annot', 			{{'triangle'}} ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			200, ...
				'n_trees', 				1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				0.05 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)

		case 8
			% B cell maturation
			% Bendall et al. 2014. “Single-Cell Trajectory Detection Uncovers Progression and Regulatory Coordination in Human B Cell Development.” Cell 157 (3). Elsevier Inc.: 714–25.
			fprintf('running treetop on B cell maturation data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'B cell'), ...
				'output_dir', 			fullfile(output_dir, 'B cell'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('PC%02d', ii), 1:10, 'unif', false)}, ...
				'extra_markers', 		{{'CD10', 'CD117', 'CD179a', 'CD179b', 'CD19', 'CD20', 'CD24', 'CD34', 'CD38', 'CD45', 'CD72', 'CD79b', ...
										'HLADR', 'IgD', 'IgM-i', 'IgM-s', 'Kappa', 'Lambda'}}, ...
				'used_cofactor', 		0, ...
				'extra_cofactor', 		5, ...
				'mat_file', 			'B cells.mat' ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size',  		1e5, ...
				'n_ref_cells',  		200, ...
				'n_trees',  			1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.2, ...
				'sigma', 				5 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 9
			% Healthy human bone marrow, scRNAseq data
			% Velten et al. 2017. “Human Haematopoietic Stem Cell Lineage Commitment Is a Continuous Process.” Nature Cell Biology 19 (4): 271–81.
			fprintf('running treetop on healthy human bone marrow, scRNAseq data (Paul et al.)\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'bone marrow scRNAseq Paul'), ...
				'output_dir', 			fullfile(output_dir, 'bone marrow scRNAseq Paul'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('diff%02d', ii), 1:8, 'unif', false);}, ...
				'used_cofactor', 		0, ...
				'mat_file', 			'paul diffusion maps.mat' ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'outlier', 				0, ...
				'threshold', 			0.2, ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			100, ...
				'n_dens_ref', 			500, ...
				'n_trees', 				1000, ...
				'layout_tree_idx', 		4, ...
				'sigma', 				1e-2 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 10
			% Healthy human bone marrow, scRNAseq data
			% Velten et al. 2017. “Human Haematopoietic Stem Cell Lineage Commitment Is a Continuous Process.” Nature Cell Biology 19 (4): 271–81.
			fprintf('running treetop on healthy human bone marrow, scRNAseq data (Velten et al.)\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'bone marrow scRNAseq Velten'), ...
				'output_dir', 			fullfile(output_dir, 'bone marrow scRNAseq Velten'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('DC%02d', ii), 1:11, 'unif', false);}, ...
				'extra_markers', 		{{'stemnet_p'}}, ...
				'mat_file', 			'velten diffusion maps, STEMNET labels.mat' ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'outlier', 				0, ...
				'threshold', 			0.1, ...
				'sample_size', 			1e5, ...
				'n_ref_cells', 			100, ...
				'n_dens_ref', 			500, ...
				'n_trees', 				1000, ...
				'layout_tree_idx', 		2, ...
				'sigma', 				1 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)
			treetop_recursive(input_struct, options_struct)

		case 11
			% Swiss roll synthetic data
			fprintf('running treetop on Swiss roll synthetic data\n')
			input_struct 	= struct( ...
				'data_dir', 			fullfile(data_dir, 'synthetic swiss roll'), ...
				'output_dir', 			fullfile(output_dir, 'synthetic swiss roll'), ...
				'used_markers', 		{arrayfun(@(ii) sprintf('S%02d', ii), 1:10, 'unif', false)}, ...
				'extra_markers', 		{{'angle'}}, ...
				'used_cofactor', 		0, ...
				'extra_cofactor', 		0, ...
				'mat_file', 			'swiss_roll.mat' ...
				);
			options_struct 	= struct( ...
				'metric_name', 			'L1', ...
				'sample_size',  		1e5, ...
				'n_ref_cells',  		200, ...
				'n_trees',  			1000, ...
				'outlier', 				0.01, ...
				'threshold', 			0.5, ...
				'sigma', 				1 ...
				);
			treetop_pre_run(input_struct, options_struct)
			treetop(input_struct, options_struct)

		otherwise
			error('invalid run_switch')

	end
end
