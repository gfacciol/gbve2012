function uu = vp(u0, m, uinit, v_fw, v_bw, params, Kf, Kb, reference_u0)
% uu = vp(u0, m, uinit, v_fw, v_bw, params , [Kf , Kb, reference_u0])
% Solves: \sum \| \nabla \partial_v u - g \| ^2  st  u|_{m^c}  = u0
% with   g = \nabla \partial_v reference_u0
%
% u0, uinit can be color images, each channels are processed independently
% if reference_u0 is not provided then it is assumed to be 0.
% Order of the dimensions in the videos matrices: [y, x, frame, channel]

algo = params.ALGO;
beta = params.BETA;
lambda = params.LAMBDA;
gamma = params.GAMMA;
interp = params.INTERP;

[nx,ny,nt,nch] = size(u0);
M = nx*ny*nt;               % number of pixels in one channel

% occlusion masks are optional 
if ~ exist('Kb','var')
    Kf = true(nx,ny,nt);
    Kb = true(nx,ny,nt);
end


%% OPERATORS: domains, derivatives, tensors, masks
% the edit domain mask
O = domain_mask_operator(m);

% masks for the even/odd frames
tmpvid = ones(nx,ny,nt);
tmpvid(:,:, 1:2:nt ) = tmpvid(:,:,[1:2:nt])*0;
O_odd  = double ( domain_mask_operator(tmpvid)  );
O_even = double ( domain_mask_operator(1-tmpvid));
clear tmpvid

% convective derivatives (fw and bw)
disp(['Using ' interp ' interpolation'])
[Jv_fw,Sf] = fw_convective_deriv( m, v_fw, interp);
[Jv_bw,Sb] = bw_convective_deriv( m, v_bw, interp);

% compose the removed convective derivatives with the occlusion masks
Kf = Kf.*Sf;
Kb = Kb.*Sb;
sKf = domain_mask_operator(Kf);  % sparse diagonal  masks
sKb = domain_mask_operator(Kb);

% build the occlusion tensor for the forward difference gradient.
k1 = cat(1, Kf([2:end],:,:), zeros(1,ny,nt)) .* Kf(:,:,:);
k2 = cat(2, Kf(:,[2:end],:), zeros(nx,1,nt)) .* Kf(:,:,:);
Kappa_f = [ domain_mask_operator(k1),  sparse(M,M) ; ...
            sparse(M,M),  domain_mask_operator(k2) ];
k1 = cat(1, Kb([2:end],:,:), zeros(1,ny,nt)) .* Kb(:,:,:);
k2 = cat(2, Kb(:,[2:end],:), zeros(nx,1,nt)) .* Kb(:,:,:);
Kappa_b = [ domain_mask_operator(k1),  sparse(M,M) ; ...
            sparse(M,M),  domain_mask_operator(k2) ];  

% spatial gradient with forward differences.
[Grad1,Grad2] = gradient_operator( m );  
Grad = [Grad1;Grad2];
clear Grad1 Grad2 Sf Sb k1 k2



%% Build the least squares system of equations
switch algo
    case {'FWBW', 'FWBW_GBC', 'BWFW', 'BWFW_GBC'}
        disp(['GBC: ' num2str(1-beta) ' fwbw + ' num2str(beta) ' bwfw'])
        A = [ sqrt(1-beta)*Kappa_f *Grad * O_even * Jv_fw; ...   % E^even GBC
              sqrt(1-beta)*Kappa_b *Grad * O_even * Jv_bw; ...
              sqrt(beta)  *Kappa_f *Grad * O_odd  * Jv_fw; ...   % E^odd  GBC
              sqrt(beta)  *Kappa_b *Grad * O_odd  * Jv_bw; ...
              sqrt(lambda)  * Grad                    ];
        disp(['lambda ' num2str(lambda) ' spatial regularization'])
        
        if gamma ~= 0                % this not part of the main algorithm
        disp(['gamma ' num2str(gamma) ' temporal regularization'])
            A2 = [sqrt(1-beta)* sKf * O_even * Jv_fw; ...    % E^even  BC
                  sqrt(1-beta)* sKb * O_even * Jv_bw; ...
                  sqrt(beta)  * sKf * O_odd  * Jv_fw; ...    % E^odd   BC
                  sqrt(beta)  * sKb * O_odd  * Jv_bw     ];
            A = [A; sqrt(gamma)*A2];       % concatenate
            clear A2
        end
       
          
    case {'FW_GBC', 'FW'}
        disp('GBC: fw')
        A = Kappa_f * Grad * Jv_fw;    % GBC
        
            
    case {'BW_GBC', 'BW'}
        disp('GBC: bw')
        A = Kappa_b * Grad * Jv_bw;    % GBC
        
        
     case {'FWBW_BC', 'BWFW_BC'}
        disp(['BC: ' num2str(1-beta) ' fwbw + ' num2str(beta) ' bwfw'])
        A = [ sqrt(1-beta)* sKf * O_even * Jv_fw; ...    % E^even  BC
              sqrt(1-beta)* sKb * O_even * Jv_bw; ...
              sqrt(beta)  * sKf * O_odd  * Jv_fw; ...    % E^odd   BC
              sqrt(beta)  * sKb * O_odd  * Jv_bw; ...
              sqrt(lambda)  * Grad                 ];       
        disp(['lambda ' num2str(lambda) ' spatial regularization'])

                                      
    case 'FW_BC'
        disp('BC: fw')
        A =    sKf * Jv_fw;      % BC
        
   
    case 'BW_BC'
        disp('BC: bw')
        A =   sKb * Jv_bw;      % BC
        
end
clear sKf O_even Jv_fw sKb  O_odd   Jv_bw Grad  Kappa_f Kappa_b

%% CALL the solver for each channel with initial condition uinit
uu=uinit;
for ch=1:nch
    % in case an attachment term is provided
    if exist('reference_u0','var')
        tmp =  reference_u0(:,:,:,ch);
        b = A * tmp(:);
    else
        b = [];
    end
    tmp = solve_least_squares_with_dirichlet_bc (A , m, u0(:,:,:,ch), ...
                                                  b, uinit(:,:,:,ch),params);
    uu(:,:,:,ch) = reshape(tmp, size(uu(:,:,:,ch)) );
end
