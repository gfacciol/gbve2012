function [Jv_fw,Sf]=fw_convective_deriv( u0, v_fw, interp)
% forward convective derivative with neumann boundaries
% u0 is needed ily for computing the size of the video
% v_fw(:,:,:,[1,2]) are the flow components 

if ~exist('interp','var')
   interp='BILINEAR';
end

[nx,ny,nt] = size(u0);
M = numel(u0);

vx_fw = v_fw(:,:,:,2); % the conventions for x,y are reversed in this file
vy_fw = v_fw(:,:,:,1);

% start preparing the forward warped image

% generate the index positions
[px,py,pt] = ind2sub( size(u0), [1:numel(u0)]);
i = sub2ind(size(u0),px,py,pt); % rows    (x,y,t)

% interpolate 
ivx = floor(vx_fw);
ivy = floor(vy_fw);
dx = vx_fw(:) - ivx(:);
dy = vy_fw(:) - ivy(:);

switch upper(interp)
case 'BICUBIC'
   s = bicubic_coeffs(dx,dy) ;
otherwise
   s = bilinear_coeffs(dx,dy) ;
end

neigh_size=size(s,1);
n1 = -floor((neigh_size-1)/2);

% compute the bilinear coefficients
%Iv = sparse(M,M);
% optimized 
ntriplets = M*neigh_size*neigh_size ;
tIv_I = zeros (ntriplets, 1) ;
tIv_J = zeros (ntriplets, 1) ;
tIv_X = zeros (ntriplets, 1) ;
current_triplet=1;

outMsk = px'*0;  % mask of points that fall outside of the image
for d1 =n1:n1+neigh_size-1
for d2 =n1:n1+neigh_size-1
    % index of the points for each stencil position
    pxivx = max(min(px'+ivx(:)+d1,nx),1);
    pyivy = max(min(py'+ivy(:)+d2,ny),1);
    pti =   max(min(pt' + 1, nt),1);
    
    % if some pixel fall outside we account for it in outMsk
    % this outmask will later be used to implement: S^_t
    outMsk = outMsk + (pxivx ~= (px'+ivx(:)+d1));
    outMsk = outMsk + (pyivy ~= (py'+ivy(:)+d2));
    outMsk = outMsk + (pti ~= (pt' + 1));
    
    % optimized sparse matrix construction
%    j = sub2ind(size(u0),pxivx,pyivy,pti); %columns (x+vx +-1,y+vy +-1,t+1)
%    Iv = Iv + sparse(i',j,s{d1+1-n1,d2+1-n1},M,M);
    j = sub2ind(size(u0),pxivx,pyivy,pti); %columns (x+vx +-1,y+vy +-1,t+1)
    tIv_I(current_triplet:current_triplet+M-1) = i';    
    tIv_J(current_triplet:current_triplet+M-1) = j;    
    tIv_X(current_triplet:current_triplet+M-1) = s{d1+1-n1,d2+1-n1};    
    current_triplet = current_triplet+M;
end
end

% OPTIMIZATION
%http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/#9
Iv = sparse( tIv_I, tIv_J, tIv_X,M,M);
clear tIv_I tIv_J tIv_X

% identity
I0 = speye(M,M);

% compute the convective derivative
Jv_fw = Iv - I0;

% remove all entries for the last frame. Going forward from if gets outside
% of the video and these cases are not handled during the construction 
% so we do it here. This is S^f_t.
tidx = find (outMsk==0);
ti = sub2ind(size(u0),px(tidx),py(tidx),pt(tidx));
Otmp = sparse(ti,ti,ones(size(ti)),M,M);
Jv_fw = Otmp * Jv_fw;

% S^f_t : the mask of removed entries
Sf = reshape((outMsk==0),size(u0));

% test 
% tmp = reshape( Jv * u0(:) , size(u0));
% %
% figure(2);
% for t=[1:nt]
%    imagesc(tmp(:,:,t))
%    colorbar
%    pause
% end
end



