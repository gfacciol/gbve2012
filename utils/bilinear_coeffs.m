function s = bilinear_coeffs(dx,dy) 
% generate bilinear coefficients for interpolating at 
% % positions (dx,dy) \in [0,1]^2

cx= [1-dx(:) , dx(:)];
cy= [1-dy(:) , dy(:)];
s={};
for x=1:2
    for y=1:2
        s{x,y}= cx(:,x) .* cy(:,y);
    end
end
