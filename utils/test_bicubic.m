function uu=testbicubic( u0, v_bw)


[nx,ny] = size(u0);
M = numel(u0);

vx_bw = v_bw(:,:,2); % the conventions for x,y are reversed in this file
vy_bw = v_bw(:,:,:,1);

% start preparing the backward warped image

% generate the index positions
[px,py] = ind2sub( size(u0), [1:M]);
i = sub2ind(size(u0),px,py); % rows    (x,y,t)

% interpolate 
ivx = floor(vx_bw);
ivy = floor(vy_bw);
dx = vx_bw(:) - ivx(:);
dy = vy_bw(:) - ivy(:);

s={};
s{1,1} = (1-dx) .* (1-dy);
s{2,1} = (dx) .* (1-dy);
s{1,2} = (1-dx) .* (dy);
s{2,2} = (dx) .* (dy);

s = bicubic_coeffs(dx,dy) ;

neigh_size=size(s,1);
n1 = -floor((neigh_size-1)/2);
n2 = n1+neigh_size-1;

% compute the bilinear coefficients
Iv = sparse(M,M);
outMsk = px'*0;  % mask of points that fall outside of the image
for d1 =n1:n2
for d2 =n1:n2
    pxivx = max(min(px'+ivx(:)+d1,nx),1);
    pyivy = max(min(py'+ivy(:)+d2,ny),1);
    
    outMsk = outMsk + (pxivx ~= (px'+ivx(:)+d1));
    outMsk = outMsk + (pyivy ~= (py'+ivy(:)+d2));
    
    j = sub2ind(size(u0),pxivx,pyivy); %columns (x+vx +-1,y+vy +-1,t-1)
    Iv = Iv + sparse(i',j,s{d1+1-n1,d2+1-n1},M,M);
end
end

uu=Iv * u0(:);
uu=reshape(uu,size(u0));