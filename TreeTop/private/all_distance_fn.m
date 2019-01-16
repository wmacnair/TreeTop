%% all_distance_fn: flexible and memory-safe distance calculation function
% expects X, Y to be formatted as n observations by d features for both
function dist = all_distance_fn(X, Y, metric)
	switch metric
		case 'L1'
			dist 	= comp_dist_L1(X, Y);

		case 'L2'
			dist 	= comp_dist_euclidean(X, Y);

		case 'angle'
			dist 	= comp_dist_angle(X, Y);

		case 'corr'
			dist 	= comp_dist_corr(X, Y);

		otherwise
			error('invalid metric selected');
	end
end

%% comp_dist_euclidean: 
function dist = comp_dist_euclidean(X, Y)
	% X has size m*d, Y n*d
	m 		= size(X, 1);
	n 		= size(Y, 1);
	dist 	= zeros(m,n);

	for ii = 1:m
		dist(ii,:) = sqrt(sum((repmat(X(ii,:),n,1) - Y).^2,2)); 
	end
end

%% comp_dist_L1: 
function dist = comp_dist_L1(X, Y)
	% set up
	m 		= size(X, 1);
	n 		= size(Y, 1);
	dist 	= zeros(m, n);

	for ii = 1:m
		dist(ii,:) = sum(abs(repmat(X(ii,:), n, 1) - Y),2); 
	end
end

%% comp_dist_corr: 
function dist = comp_dist_corr(X, Y)
	corr 	= calc_corr(X,Y);
	dist 	= 1-corr; 
end

%% comp_dist_angle: 
function dist = comp_dist_angle(X, Y)
	% L2-normalization
	X 	= bsxfun(@times, X, 1./sqrt(sum(X.^2, 2)));
	Y 	= bsxfun(@times, Y, 1./sqrt(sum(Y.^2, 2)));

	% dot_prod 	= X*Y';
	% dist 		= acos(dot_prod)/pi;

	% set up
	m 		= size(X, 1);
	n 		= size(Y, 1);
	dist 	= zeros(m, n);

	for ii = 1:m
		dot_prod 		= X(ii, :) * Y';
		dist(ii, :) 	= real(acos(dot_prod)/pi);
	end
end

%% comp_dist_abs_corr: 
function dist = comp_dist_abs_corr(X, Y)
	corr 		= calc_corr(X,Y);
	dist 		= 1-abs(corr); 
end

%% calc_corr: 
function [corr] = calc_corr(X, Y)
	% set up
	m 		= size(X, 1);
	n 		= size(Y, 1);
	corr 	= zeros(m, n);

	% zero-mean
	X 		= bsxfun(@minus, X, mean(X, 1));
	Y 		= bsxfun(@minus, Y, mean(Y, 1));

	% L2-normalization
	X 		= bsxfun(@times, X, 1./sqrt(sum(X.^2, 1)));
	Y 		= bsxfun(@times, Y, 1./sqrt(sum(Y.^2, 1)));

	% calculation correlation
	for ii = 1:m
		corr(ii,:) 	= X(ii,:) * Y';
	end
end
