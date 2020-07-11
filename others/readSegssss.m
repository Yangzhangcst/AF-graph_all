function [segs] = readSegssss(bsdsRoot,present,iid)

files = fullfile(bsdsRoot,'groundTruth',[iid,'.mat']);
segs_gt = load(files);
segs = cell(1,length(segs_gt.groundTruth));
for i = 1:length(segs_gt.groundTruth)
    segs{1,i} = double(segs_gt.groundTruth{1,i}.Segmentation);
end
% segs = cell(1,length(segs_gt));
% for i = 1:length(segs_gt)
%     segs{1,i} = double(segs_gt.newseg);
% end
if 0
    files = fullfile(bsdsRoot,'groundTruth_layers',[iid,'.mat']);
    segs_gt = load(files);
    segs = cell(1,length(segs_gt.groundTruth));
    for i = 1:length(segs_gt.groundTruth)
        segs{1,i} = double(segs_gt.groundTruth{1,i}.Segmentation);
    end
end
