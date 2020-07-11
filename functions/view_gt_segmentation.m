function [gt_imgs gt_cnt] = view_gt_segmentation(bsdsRoot,img,present,out_path,only_name,para)

if para.Nimgs == 300
    gt_imgs = readSegs(bsdsRoot,'color',str2num(only_name));% BSDS300
elseif para.Nimgs == 500
    gt_imgs = readSegss(bsdsRoot,present,only_name);% BSDS500
elseif para.Nimgs == 591
    gt_imgs = readSegsss(bsdsRoot,present,only_name);% MSRC
else
    gt_imgs = readSegssss(bsdsRoot,present,only_name); % Others
end
    
gt_path = fullfile(out_path, 'gt'); 
if  ~exist(gt_path), mkdir(gt_path); end

[X Y Z] = size(img); 

gt_cnt = [];
for i=1:size(gt_imgs,2)
    if para.save == 1
        [imgMasks,segOutline,imgMarkup]=segoutput(img,gt_imgs{i});
        imwrite(segOutline,fullfile(gt_path, [only_name, '_outline_', int2str(i), '.bmp'])); 
        imwrite(imgMarkup,fullfile(gt_path, [only_name, '_', int2str(i), '.bmp'])); 
        disp_img = zeros(X,Y,3); L = max(gt_imgs{i}); 
        for k=1:L
            idx = find(imgMasks==k);
            for j=1:3
                tmp = disp_img(:,:,j); tmp1 = img(:,:,j);
                tmp(idx) = mean(tmp1(idx)); 
                disp_img(:,:,j) = tmp; 
                clear tmp tmp1;
            end
            clear idx;
        end
        [~,~,imgMasks]=segoutput(disp_img,gt_imgs{i});
        imwrite(imgMasks,fullfile(gt_path, [only_name, '_mask_', int2str(i), '.bmp'])); 
    end
    if iscell(gt_imgs)
        gt_cnt(i) = max(gt_imgs{i}(:));
    else
        gt_cnt(i) = max(gt_imgs(:));
    end
end
