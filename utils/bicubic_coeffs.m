function s = bicubic_coeffs(dx,dy) 
% generate bicubic coefficients for interpolating at 
% positions (dx,dy) \in [0,1]^2
%
% s{1,1} s{2,1} s{3,1}    s{4,1}
%  o      o      o      o 
%  
% 
%                         s{4,2} 
%  o      X--dx  o      o
%         |
%         dy  x 
%                         s{4,3} 
%  o      o      o      o
% 
%  
%                         s{4,4} 
%  o      o      o      o
%  

cx= keys(dx(:));
cy= keys(dy(:));

s={};
for x=1:4
    for y=1:4
        s{x,y}= cx(:,x) .* cy(:,y);
    end
end
end

function c = keys (t,a)
% c = keys(t,[a=-0.5])
% coefficients for 1D cubic interpolant (Keys's) are stored as stencil
% the center sample is at c[2] and the sample corresponding
% to the left sample of x (i.e. x-1) is c[1]
% c[1]   c[2]    c[3]    c[4]
%        t=0 <--> t=1
% x-1     x       x+1

if ~exist('a','var')
    a=-0.5;
end

t  = t(:);
t2 = t(:) .* t(:);
at = a * t(:);
c(:,4) = a * t2 .* (1.0 - t);
c(:,3) = (2.0 * a + 3.0 - (a + 2.0) .* t) .* t2 - at;
c(:,2) = ((a + 2.0) .* t - a - 3.0) .* t2 + 1.0;
c(:,1) = a * (t - 2.0) .* t2 + at;
end

