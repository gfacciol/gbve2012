%%% PARAMETERS
clear
addpath('utils')
addpath('utils/multi')

FRAMENUM    =  [2:60]; %[1,2,54,55]
IMAGES_PATT = 'data/mesa/frames/Images_%03d.png';
MASKS_PATT  = 'data/mesa/mask/Mask_%03d.png';
FLOW_FW_PATT= 'data/mesa/OFD/HOMOGRfw%03d.mat';
FLOW_BW_PATT= 'data/mesa/OFD/HOMOGRbw%03d.mat';
%FLOW_FW_PATT= 'data/mesa/OFD/SUN_fw%03d.mat';
%FLOW_BW_PATT= 'data/mesa/OFD/SUN_bw%03d.mat';
%OCC_FW_PATT = 'data/cloth/OFD/occ_gbve_fw%03d.png'; 
%OCC_BW_PATT = 'data/cloth/OFD/occ_gbve_bw%03d.png'; 
OUT_DIR     = 'data/mesa/results/';
LIDS        = 'FIRST';  % 'FIRST' 'LAST', 'BOTH'

params.ALGO='FWBW_BC';        %'FWBW_GBC' 'FW_GBC' 'BW_GBC' 'FWBW_BC' 'FW_BC' 'BW_BC'
                        % GBC=Global Brightness Change model
                        % BC =Brightness Constancy model

params.BETA=1;              % FBBW/BWFW mix    ONLY FOR FWBW_*
params.LAMBDA=0.00;            % SPATIAL REGULARITY 0.02    ONLY FOR FWBW_* 
params.GAMMA =0.00;            % TEMPORAL REGULARITY USUALLY 0, ONLY for FWBW_GBC 
params.INTERP='BICUBIC';      % INTERPOLATION: 'BILINEAR' 'BICUBIC'
params.CGprec=0.001;

%%% PRE ALLOCATE VIDEOS
%  dimensions: [y, x, frame, channel]
[nx,ny,nch]    = size(imread(sprintf(IMAGES_PATT , FRAMENUM(1) )));
%u0    = zeros( nx, ny, numel(FRAMENUM), nch);
u0    = zeros( nx, ny, numel(FRAMENUM), 1);
m     = false( nx, ny, numel(FRAMENUM), 1);
v_fw  = zeros( nx, ny, numel(FRAMENUM), 2);
v_bw  = zeros( nx, ny, numel(FRAMENUM), 2);
occ_fw= true ( nx, ny, numel(FRAMENUM), 1);
occ_bw= true ( nx, ny, numel(FRAMENUM), 1);

%%% LOAD VIDEOS AND FLOWS
j=1;
for i =  FRAMENUM 
    disp (['loading frame ', num2str(i), '...'])
    % frames
    Im = sprintf(IMAGES_PATT , i);
    Im = imread(Im);
    if(size(Im,3)==3)            % FOR THE MOMENT ONLY BW
        %u0(:,:,j,:)=double((Im)); 
        u0(:,:,j,:)=double(rgb2gray(Im)); 
    else
        u0(:,:,j)=double(Im(:,:,1)); 
    end   
    
    % masks
    Im = sprintf(MASKS_PATT , i);
    Im = imread(Im);
    m(:,:,j) = logical(Im(:,:,1));
    
    % flows FW
    Im = sprintf(FLOW_FW_PATT , i);
    if(exist(Im,'file') )
        load(Im)
    else
        uv = repmat ( m(:,:,j)*0,[1,1,2]);
        disp(['not found (setting to 0):', Im])
    end
    v_fw(:,:,j,:) = uv(:,:,:);
    
    % flows BW
    Im = sprintf(FLOW_BW_PATT , i);
    if(exist(Im,'file') )
        load(Im)
    else
        uv = repmat ( m(:,:,j)*0,[1,1,2]);
        disp(['not found (setting to 0):', Im])
    end
    v_bw(:,:,j,:) = uv(:,:,:);    
    
%     % occlusions FW BW
%     Im = imread(sprintf(OCC_FW_PATT , i));
%     occ_fw(:,:,j) = imerode(double(Im(:,:,1)) ~= 255,ones(1)); 
%     
%     Im = imread(sprintf(OCC_BW_PATT , i));
%     occ_bw(:,:,j) = imerode(double(Im(:,:,1)) ~= 255,ones(1));    
    j=j+1;
end

[nx,ny,nt,nch] = size(u0);
M = numel(u0);



%%% crop interest region +  something extra (10 pixels)
cor   = compute_ROI( m , 10) ;

sstep = 1;
u0    =    u0(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:,:);
m     =     m(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:);
v_fw  =  v_fw(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:,:)/sstep;
v_bw  =  v_bw(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:,:)/sstep;
occ_fw=imerode(occ_fw,ones(sstep));
occ_bw=imerode(occ_bw,ones(sstep));
occ_fw=occ_fw(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:);
occ_bw=occ_bw(cor(1):sstep:cor(3), cor(2):sstep:cor(4),:);



%%% set mask at the lid
switch LIDS
    case 'FIRST'
        m(:,:,1) = 0;
    case 'LAST'
        m(:,:,nt) = 0;        
    case 'BOTH'
        m(:,:,1) = 0;
        m(:,:,nt) = 0;        
end





%% show input
clf; colormap gray;
for t= 1:nt
    subplot(1,3,1);  imagesc( squeeze(u0(:,:,t,:))/255, [0,1])
    subplot(1,3,2);  imagesc( m(:,:,t))
    subplot(1,3,3);  imagesc( occ_fw(:,:,t))
    drawnow;
end


%%
% tic
% uu=u0;
% for t=1:2:nt-1
%     mv=m;
%     mv(:,:,t) = 0;
%     tmp = vp(uu(:,:,t:t+2,:), mv(:,:,t:t+2), uu(:,:,t:t+2,:),...
%         v_fw(:,:,t:t+2,:), v_bw(:,:,t:t+2,:), ALGO, BETA, LAMBDA, GAMMA );
%     uu(:,:,t:t+2,:) =  reshape(tmp,size(uu(:,:,t:t+2,:)));
% end
% toc
% u0=uu;

%%% call the method
tic
uu = vp(u0, m, u0, v_fw, v_bw, params , occ_fw, occ_bw);
toc


%% show result
clf; colormap gray;
mkdir(OUT_DIR);    
uuf = filterHF(uu);
for t= 1:nt
    subplot(1,3,1);  imagesc( squeeze(u0(:,:,t,:))/255, [0,1])
    subplot(1,3,2);  imagesc( min(max(squeeze(uu(:,:,t,:))/255,0),1), [0,1])
    subplot(1,3,3);  imagesc( min(max(squeeze(uuf(:,:,t,:))/255,0),1), [0,1])
    drawnow;    pause(0.1);

%    imwrite(uint8 ( min(max(squeeze(uu(:,:,t,:)),0),255)), sprintf('%s/out%s_b%g_l%g_g%g_%03d.png', OUT_DIR,params.ALGO, params.BETA, params.LAMBDA, params.GAMMA, t) )
%    imwrite(uint8 ( min(max(squeeze(uu(:,:,t,:)),0),255)), sprintf('%s/out_filt_%s_b%g_l%g_g%g_%03d.png', OUT_DIR,params.ALGO, params.BETA, params.LAMBDA, params.GAMMA, t) )
end
