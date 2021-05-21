function [Jv_bw,Sb]=bw_convective_deriv( u0, v_bw)
% backward convective derivative with neumann boundaries
% u0 is needed ily for computing the size of the video
% v_bw(:,:,:,[1,2]) are the flow components 

[nx,ny,nt] = size(u0);
M = numel(u0);

vx_bw = v_bw(:,:,:,2); % the conventions for x,y are reversed in this file
vy_bw = v_bw(:,:,:,1);

% start preparing the backward warped image

% generate the index positions
[px,py,pt] = ind2sub( size(u0), [1:M]);
i = sub2ind(size(u0),px,py,pt); % rows    (x,y,t)

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

% compute the bilinear coefficients
Iv = sparse(M,M);
outMsk = px'*0;  % mask of points that fall outside of the image
for d1 =0:1
for d2 =0:1
    pxivx = max(min(px'+ivx(:)+d1,nx),1);
    pyivy = max(min(py'+ivy(:)+d2,ny),1);
    pti =   max(min(pt' - 1, nt),1);
    
    outMsk = outMsk + (pxivx ~= (px'+ivx(:)+d1));
    outMsk = outMsk + (pyivy ~= (py'+ivy(:)+d2));
    outMsk = outMsk + (pti ~= (pt' - 1));
    
    j = sub2ind(size(u0),pxivx,pyivy,pti); %columns (x+vx +-1,y+vy +-1,t-1)
    Iv = Iv + sparse(i',j,s{d1+1,d2+1},M,M);
end
end

% identity
I0 = speye(M,M);

% compute the convective derivative
Jv_bw = I0 - Iv;

% remove all entries for the last frame. Going forward from if gets outside
% of the video and these cases are not handled during the construction 
% so we do it here. This is S^b_t.
tidx = find (outMsk==0);
ti = sub2ind(size(u0),px(tidx),py(tidx),pt(tidx));
Otmp = sparse(ti,ti,ones(size(ti)),M,M);
Jv_bw = Otmp * Jv_bw;

% S^b_t : the mask of removed entries
Sb = reshape((outMsk==0),size(u0));

% test 
% tmp = reshape( Jv * u0(:) , size(u0));
% %
% figure(2);
% for t=[1:nt]
%    imagesc(tmp(:,:,t))
%    colorbar
%    pause
% end

