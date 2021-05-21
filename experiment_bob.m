%%% PARAMETERS
clear
addpath('utils')

FRAMENUM    =  [2:38]; %[1,2,54,55]
IMAGES_PATT = 'data/bob/input_color/frame%03d.png';
MASKS_PATT  = 'data/bob/masks/mask_%03d.png';
FLOW_FW_PATT= 'data/bob/OFD/fw_%03d.mat';
FLOW_BW_PATT= 'data/bob/OFD/bw_%03d.mat';
%OCC_FW_PATT = 'data/bob/OFD/occ_gbve_fw%03d.png'; 
%OCC_BW_PATT = 'data/bob/OFD/occ_gbve_bw%03d.png'; 
OUT_DIR     = 'data/bob/results/';
LIDS        = 'FIRST';  % 'FIRST' 'LAST', 'BOTH'

ALGO='FWBW_GBC';        %'FWBW_GBC' 'FW_GBC' 'BW_GBC' 'FWBW_BC' 'FW_BC' 'BW_BC'
                        % GBC=Global Brightness Change model
                        % BC =Brightness Constancy model

BETA=0.95;              % FBBW/BWFW mix    ONLY FOR FWBW_*
LAMBDA=0.00;            % SPATIAL REGULARITY 0.02    ONLY FOR FWBW_* 
GAMMA =0.00;            % TEMPORAL REGULARITY USUALLY 0, ONLY for FWBW_GBC 

%%% PRE ALLOCATE VIDEOS
%  dimensions: [y, x, frame, channel]
[nx,ny,nch]    = size(imread(sprintf(IMAGES_PATT , FRAMENUM(1) )));
u0    = zeros( nx, ny, numel(FRAMENUM), nch);
%u0    = zeros( nx, ny, numel(FRAMENUM), 1);
m     = zeros( nx, ny, numel(FRAMENUM), 1);
v_fw  = zeros( nx, ny, numel(FRAMENUM), 2);
v_bw  = zeros( nx, ny, numel(FRAMENUM), 2);
occ_fw= ones ( nx, ny, numel(FRAMENUM), 1);
occ_bw= ones ( nx, ny, numel(FRAMENUM), 1);

%%% LOAD VIDEOS AND FLOWS
j=1;
for i =  FRAMENUM 
    disp (['loading frame ', num2str(i), '...'])
    % frames
    Im = sprintf(IMAGES_PATT , i);
    Im = imread(Im);
    if(size(Im,3)==3)            % FOR THE MOMENT ONLY BW
        u0(:,:,j,:)=double((Im)); 
        %u0(:,:,j,:)=double(rgb2gray(Im)); 
    else
        u0(:,:,j)=double(Im(:,:,1)); 
    end   
    
    % masks
    Im = sprintf(MASKS_PATT , i);
    Im = imread(Im);
    m(:,:,j) = double(Im(:,:,1));
    
    % flows FW
    Im = sprintf(FLOW_FW_PATT , i);
    if(exist(Im,'file') )
        load(Im)
    else
        uv = repmat ( m(:,:,j)*0,[1,1,2]);
    end
    v_fw(:,:,j,:) = uv(:,:,:);
    
    % flows BW
    Im = sprintf(FLOW_BW_PATT , i);
    if(exist(Im,'file') )
        load(Im)
    else
        uv = repmat ( m(:,:,j)*0,[1,1,2]);
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
figure(1);    colormap gray;
for t= 1:nt
    subplot(1,3,1);  imagesc( squeeze(u0(:,:,t,:))/255, [0,1])
    subplot(1,3,2);  imagesc( m(:,:,t))
    subplot(1,3,3);  imagesc( occ_fw(:,:,t))
    drawnow;
end


%%% call the method
tic
uu = vp(u0, m, u0, v_fw, v_bw, ALGO, BETA, LAMBDA, GAMMA, occ_fw, occ_bw);
toc


%% show result
figure(1);    colormap gray;
mkdir(OUT_DIR);    
uuf = filterHF(uu);
for t= 1:nt
    subplot(1,3,1);  imagesc( squeeze(u0(:,:,t,:))/255, [0,1])
    subplot(1,3,2);  imagesc( min(max(squeeze(uu(:,:,t,:))/255,0),1), [0,1])
    subplot(1,3,3);  imagesc( min(max(squeeze(uuf(:,:,t,:))/255,0),1), [0,1])
    drawnow;    pause(0.1);

    imwrite(uint8 ( min(max(squeeze(uu(:,:,t,:)),0),255)), sprintf('%s/out%s_b%g_l%g_g%g_%03d.png', OUT_DIR,ALGO, BETA, LAMBDA, GAMMA, t) )
    imwrite(uint8 ( min(max(squeeze(uu(:,:,t,:)),0),255)), sprintf('%s/out_filt_%s_b%g_l%g_g%g_%03d.png', OUT_DIR,ALGO, BETA, LAMBDA, GAMMA, t) )
end
