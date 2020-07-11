clc; clear; close all;

%% set path
addpath 'others'
addpath 'evals'
addpath 'KSC'
addpath 'SSC'
addpath 'functions'
addpath 'msseg'
addpath 'APC'
addpath 'Graph_based_segment'


%% set parameters for bipartite graph
para.alpha = 0.001; % affinity between pixels and superpixels
para.beta  =  20;   % scale factor in superpixel affinity
para.nb = 1; % number of neighbors for superpixels
para.rho = 1;
para.alphak = 1e-5;
para.betak = 1e-6;
para.Nimgs = 591; % number of images in MSRC
para.save = 0;

%% read numbers of segments used in the paper 
bsdsFile = 'MSRC';
bsdsRoot = ['./database/',bsdsFile];
saveRoot = 'results';
fid = fopen(fullfile('results_searchN',bsdsFile,'Nsegs.txt'),'r');
line = 1;
while feof(fid) == 0    
    BSDS_INFO{line,1} = deblank(fgetl(fid)); 
    line = line+1;
end
fclose(fid);

%% PRI,VoI,GCE,BDE.
PRI_all = zeros(para.Nimgs,1);
VoI_all = zeros(para.Nimgs,1);
GCE_all = zeros(para.Nimgs,1);
BDE_all = zeros(para.Nimgs,1);

%% Settings
% setup
nGCluster = 3; % number of subjects.
% dimension reduction
reduceDimension = @(data) dimReduction_PCA(data, 0);
% normalization
normalizeColumn = @(data) cnormalize_inplace(data);
% representation
buildRepresentation = @(data) SP_mat_func(data, 3, 1e-6); % second parameter is sparsity
% spectral clustering   
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');

%%
Nseg_save = [];
for idxI = 1:para.Nimgs
    
   % read number of segments
    S = regexp(BSDS_INFO{idxI},'\s+','split');tic;
    if length(S)>2 
        img_name = S{2};Nseg = str2double(S{3});
    else
        img_name = S{1};Nseg = str2double(S{2});
    end
    img_loc = fullfile(bsdsRoot,'images',[img_name,'.bmp']);
    present = 'image';
    img = imread(img_loc); [X,Y,~] = size(img);
    out_path = fullfile(saveRoot,bsdsFile,img_name);
    if ~exist(out_path,'dir'), mkdir(out_path); end
%     
    % generate superpixels
    [para_MS, para_FH] = set_parameters_oversegmentation(img);
    [seg,labels_img,seg_vals,seg_lab_vals,seg_edges,seg_img] = make_superpixels(img,para_MS,para_FH);
    
    %% construct graph
        Np = X*Y;   Nsp = 0;
    for k = 1:length(seg)
        Nsp = Nsp + size(seg{k},2);
    end

    W_Y = sparse(Nsp,Nsp);
    edgesXY = [];
    j = 1;
    for k = 1:length(seg) % for each over-segmentation
        % for each over-segmentation
        feature = seg_lab_vals{k};
        % subspace-preserving representation
        tmp1 = reduceDimension(feature);
        tmp1 = normalizeColumn(tmp1);
        R = buildRepresentation(tmp1');
        R(1:length(feature)+1:end) = 0;
        A = abs(R) + abs(R)';
        nGCluster(idxI,k) = APclustering(feature);
        index_tmp = genLabel(A, nGCluster(idxI,k));

        % superpixel division
        local_nodes  = find(index_tmp == mode(index_tmp));
        global_nodes = find(index_tmp ~= mode(index_tmp));

        feature(:,all(feature == 0, 1))=[];
        [fm,fn] = size(feature);
        feature=(feature-repmat(mean(feature),fm,1))./repmat(std(feature),fm,1); 

        % adjacency-graph over all nodes
        w = makeweights(seg_edges{k},feature,para.beta);
        W_local = adjacency(seg_edges{k},w);
        W = W_local;

        % KSC-graph over global_nodes
        W_KSC = KSCGRAPH(feature,para); 
        W = assignGraphValue(W,W_KSC,global_nodes);
        W = sparse(W);

        Nk = size(seg{k},2); % number of superpixels in over-segmentation k
        W_Y(j:j+Nk-1,j:j+Nk-1) = prune_knn(W,para.nb);

        % affinities between pixels and superpixels
        for i = 1:Nk
            idxp = seg{k}{i}; % pixel indices in superpixel i
            Nki = length(idxp);
            idxsp = j + zeros(Nki,1);
            edgesXY = [edgesXY; [idxp, idxsp]];
            j = j + 1;
        end
    end
    W_XY = sparse(edgesXY(:,1),edgesXY(:,2),para.alpha,Np,Nsp);
    % affinity between a superpixel and itself is set to be the maximum 1.
    W_Y(1:Nsp+1:end) = 1;  B = [W_XY;W_Y];   
    
    % save over-segmentations
    view_oversegmentation(labels_img,seg_img,out_path,img_name,para.save);    
    
    % Transfer Cut
    label_img = Tcut(B,Nseg,[X,Y]); clear B; ti = toc;

    % save segmentation
    view_segmentation(im2double(img),label_img(:),out_path,img_name,para.save);
    
    % evaluate segmentation
    [gt_imgs, gt_cnt] = view_gt_segmentation(bsdsRoot,im2double(img),present,out_path,img_name,para);
    out_vals = eval_segmentation(label_img,gt_imgs); clear label_img gt_imgs;
    fprintf('%6s: %2d %9.6f, %9.6f, %9.6f, %9.6f %.4fs\n', img_name, Nseg,...
        out_vals.PRI, out_vals.VoI, out_vals.GCE, out_vals.BDE, ti);
    PRI_all(idxI) = out_vals.PRI;
    VoI_all(idxI) = out_vals.VoI;
    GCE_all(idxI) = out_vals.GCE;
    BDE_all(idxI) = out_vals.BDE;
end
fprintf('Mean: %14.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));

fid_out = fopen(fullfile(saveRoot,bsdsFile,'evaluation.txt'),'w');
for idxI=1:para.Nimgs
    S = regexp(BSDS_INFO{idxI},'\s+','split');
    if length(S)>2 
        img_name = S{2};Nseg = str2double(S{3});
    else
        img_name = S{1};Nseg = str2double(S{2});
    end
    fprintf(fid_out,'%6s %9.6f, %9.6f, %9.6f, %9.6f \n', img_name,...
        PRI_all(idxI), VoI_all(idxI), GCE_all(idxI), BDE_all(idxI));
end
fprintf(fid_out,'Mean: %10.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all),...
    mean(VoI_all), mean(GCE_all), mean(BDE_all));
fclose(fid_out);

