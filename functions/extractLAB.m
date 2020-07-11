function     LAB_features = extractLAB(seg_lab_vals,labels_img)
for m =1:length(labels_img)
    [labels_num(m) ,~ ] = max(max(labels_img{m})); %Ѱ�ҵ�i����ǩ���������
     img_seg_lab= cell2mat(seg_lab_vals(m));
     A_label = isnan(img_seg_lab);
      img_seg_lab(A_label==1) =0;
    for i = 1: labels_num(m)  % ���η���ÿ������        
        LAB_features1(m,i)= {img_seg_lab(i,:)};
    end
end

P_layNum = size(LAB_features1,1);%�������
for j = 1: P_layNum
    labFeatures_lay = LAB_features1(j,:);%�õ�һ�����ʵ��
    labFeatures_lay(cellfun(@isempty,labFeatures_lay))=[];%ȥ����Ԫ��
    S_layNum = size(labFeatures_lay,2);%����������
    TrainData_orig = cell2mat(labFeatures_lay);%����
    TrainData = reshape(TrainData_orig,3,S_layNum);%������������
%     for i = 1: S_layNum %�������У���λ��
%         TrainData(:,i) = TrainData(:,i)/(sqrt(sum(TrainData(:,i).^2))+eps);   
%     end

    feature = TrainData';
    feature(:,all(feature == 0, 1))=[];
    [fm,fn] = size(feature);
    feature=(feature-repmat(mean(feature),fm,1))./repmat(std(feature),fm,1); 
    LAB_features{1,j} = feature;
end