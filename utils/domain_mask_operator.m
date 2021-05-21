function  O = domain_mask_operator( mask )
% generate diagonal domain mask operator from the mask in the input
M = numel(mask);
O = spdiags( mask(:) > 0 ,0,M,M);

% OLD AND SLOW
%
%[px,py,pt] = ind2sub( size(mask), [1:numel(mask)]);
%i = sub2ind(size(mask),px,py,pt); % rows    (x,y,t)
%M = numel(mask);
%O = sparse(i,i,(mask(:) > 0),M,M);
