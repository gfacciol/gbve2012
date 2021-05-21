function [K1,K2]=gradient_operator( u0 )
% returns the vertical and horizontal derivative operators with neumann
% boundary conditions for images of size(u0)

[nx,ny,nt] = size(u0);
M = numel(u0);

%%% start preparing the gradient operator ( neumann at the boundaries ) 
% generate the index positions
[px,py,pt] = ind2sub( size(u0), [1:numel(u0)]);
i = sub2ind(size(u0),px,py,pt); % rows    (x,y,t)
j = sub2ind(size(u0),min(px+1,nx) ,py,pt); % columns  (x+1,y,t) 
k = sub2ind(size(u0),px,min(py+1,ny) ,pt); % columns  (x,y+1,t) 
s = ones(numel(i),1);

I0 = speye(M,M);

tmp = sparse(i,j,s,M,M);
K1 =  tmp - I0;

tmp = sparse(i,k,s,M,M);
K2 =  tmp - I0;

% the spatial gradient operator
K = [K1;K2];

% test it 
%t = reshape( K * u0(:) , [nx,ny,nt*2]);
