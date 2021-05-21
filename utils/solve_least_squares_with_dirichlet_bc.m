function xx = solve_least_squares_with_dirichlet_bc (A , msk, x0, b, xinit, params)
% xx = solve_linear_system_with_dirichlet_bc (A , msk, x0, b, xinit, [params])
% solves min | A x -  b |^2  , st x | msk^c = x0
% msk is 0 for the variables we want to keep fixed
% xinit is the initial condition for x

% Since x0 is fixed we perform the change of variables
% x =  M y +  x0, where M is the restriction mask 
% and rewrite the equivalent unconstrained system: 
% min_y | A M y + A x0 -  b |^2
% whose minumum satisfies 
% (M^T  A^T  A  M)  y  =  - M^T  A^T  (A x0 - b)
% the solution of the initial problem is recovered as:  x = y + x0


if ~ exist('b','var') || numel(b) == 0
    b = zeros(size(A,1),1);
end

if ~ exist('xinit','var')
    xinit = x0(:);
end

% generate a diagonal operator from the mask
M = domain_mask_operator( msk(:) );

% optimized 
AMt = (A*M)';
b1 =   AMt * ( - A * x0(:)  + b(:));
A1 =   AMt * AMt';
clear AMt

% OLD AND SLOW
%b1 =  -M' * (A' * A) * x0(:)  + M' * A' * b(:);
%A1 =   M' * (A' * A) * M ;


% the variable change also applies to the initial condition 
xx = xinit(:)-x0(:);
if exist('params','var') && isfield(params,'CGprec') 
    xx= pcg(A1,b1,params.CGprec,10000,[],[],xx);
else
    xx= pcg(A1,b1,0.000000001,10000,[],[],xx);
end
% xx= pcg(A1-0.2*speye(size(A1)),b1,0.00000001,10000,[],[],xx);
% for r=[0.2, 0.1, 0.05,0.02,0.01, 0.001, 0.0001, 0.00001]
%     xx= pcg(A1-r*speye(size(A1)),b1,0.00000001,10000,[],[],xx(:));
% end
xx = M * xx + x0(:);

