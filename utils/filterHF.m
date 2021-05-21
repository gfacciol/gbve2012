function out = filterHF(in)
% Hipomode filter of a video % G.F.

[nx,ny,nt,nch] = size(in);
% construct the filter
t=zeros(2*nx,2*ny);
t(1,1)=1/4;
t(2*nx,1)=1/4;
t(1,2*ny)=1/4;
t(2*nx,2*ny)=1/4;
filter = abs(fft2(t));
filter = filter(1:nx,1:ny);

% apply with dct to handle periodization
out=in;
for i=1:nt
    for j=1:nch
        out(:,:,i,j)= idct2(dct2(in(:,:,i,j)).*filter);
    end
end
