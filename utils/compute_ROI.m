function  corners = compute_ROI( mask , band )
% corners = compute_ROI( mask , [ band ] )
% mask is a binary mask
% band is the with of a security band
% returns [minx, miny, maxx, maxy]

if ~exist('band','var')
    band=10;
end

[nx,ny,nt] = size(mask);
M = numel(mask);

[px,py,pt] = ind2sub( size(mask), find(mask>0) );
dx = min(max(px(:))+band,nx);
cx = max(min(px(:))-band,1);
dy = min(max(py(:))+band,ny);
cy = max(min(py(:))-band,1);
 
corners= [cx, cy, dx, dy];
