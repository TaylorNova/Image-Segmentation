function [cluster,border] = mySLIC(ima,clusternum)
%%
%Ԥ����
I = ima;%ԭͼ
I = im2double(I);
I_lab = rgb2lab(I);%LABͼ��
[height,width,channel] = size(ima);
mvalue = 20;%mֵ
S = round((height*width/clusternum)^0.5);%Sֵ�ļ���
%Gabor�˲�
wavelength = 2.^(3:5) * 3;
orientation = 0:90:90;
g = gabor(wavelength,orientation); %2 * 3 = 6��
I = rgb2gray(im2single(ima));
gabormag = imgaborfilt(I,g);
for i = 1:6
    maxvalue = max(max(gabormag(:,:,i)));
    minvalue = min(min(gabormag(:,:,i)));
    len = maxvalue - minvalue;
    gabormag(:,:,i) = (gabormag(:,:,i) - minvalue)/len;%��һ��
end
%%
%�������ĵ��ʼ��
center_points = create_points(height,width,clusternum);
%��ʼ���ݶȾ�ƫ
I_grad = gradient(I_lab);%�ݶ�
I_grad_abs = abs(I_grad);%�ݶȷ���
for i = 1:clusternum
    x = center_points(i,1);
    y = center_points(i,2);
    roundpoints = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y;x,y+1;x+1,y-1;x+1,y;x+1,y+1];%3*3��Χ�ڵĵ�
    mingrad = 10000000;%��Χ����С�ݶ�
    for k = 1:9
        c_point = correct_points(roundpoints(k,:),height,width);
        if I_grad_abs(c_point(1),c_point(2)) < mingrad
            mingrad = I_grad_abs(c_point(1),c_point(2));
        end
    end
    for k = 1:9 %��С�ݶȶ�Ӧ��
        c_point = correct_points(roundpoints(k,:),height,width);
        if I_grad_abs(c_point(1),c_point(2)) == mingrad
            break;
        end
    end
    c_point = correct_points(roundpoints(k,:),height,width);
    center_points(i,1) = c_point(1);
    center_points(i,2) = c_point(2);
end
%%
%ѭ������
last_center_points = center_points;%��һ�εľ������ĵ�
cluster_array = -1*ones(height,width);%�����ǩ��ʼ��Ϊ-1
distance_array = inf*ones(height,width);%�����ʼ��Ϊ����
while true
    %������룬���¾����ǩ
    for i = 1:clusternum
        %x��y��������Χ
        x_range = [max(center_points(i,1) - S,1),min(center_points(i,1) + S,height)];
        y_range = [max(center_points(i,2) - S,1),min(center_points(i,2) + S,width)];
        for m = x_range(1):x_range(2)
            for n = y_range(1):y_range(2)
                centerlab = I_lab(center_points(i,1),center_points(i,2),:);%���ĵ��labֵ
                templab = I_lab(m,n,:);%��ǰ���labֵ
                lab_dis = (centerlab(1) - templab(1))^2 + (centerlab(2) - templab(2))^2 + (centerlab(3) - templab(3))^2;%lab�����ƽ��
                space_dis = (center_points(i,1) - m)^2 + (center_points(i,2) - n)^2;%�ռ�����ƽ��
                sumdis = 0;
                for v = 1:6
                    sumdis = sumdis + (gabormag(m,n,v) - gabormag(center_points(i,1),center_points(i,2),v))^2;
                end
                gabor_dis = sumdis/6;%gabor�˲������ƽ��
                total_dis = (lab_dis/(mvalue^2) + space_dis/(S^2) + gabor_dis)^0.5;%�ܾ���
                if total_dis < distance_array(m,n)
                    distance_array(m,n) = total_dis;
                    cluster_array(m,n) = i;
                end
            end
        end
    end
   %���¾������ĵ�
   count_sum = zeros(clusternum,2); %�����x��y�����ܺ�
   countnum = zeros(clusternum,1); %���������
   for i = 1:clusternum
        for m = 1:height
            for n = 1:width
                if cluster_array(m,n) == i
                    count_sum(i,1) = count_sum(i,1) + m;
                    count_sum(i,2) = count_sum(i,2) + n;
                    countnum(i) =  countnum(i) + 1;
                end
            end
        end
   end
   %��ֵ�����¾������ĵ�
   for i = 1:clusternum
        center_points(i,1) = round(count_sum(i,1)/countnum(i));
        center_points(i,2) = round(count_sum(i,2)/countnum(i));
   end
   %����в�ж��Ƿ����ѭ��
   dif = center_points - last_center_points;
   sumdif = 0;
   for i = 1:clusternum
       sumdif = sumdif + (dif(i,1)^2 + dif(i,2)^2)^0.5;
   end
   res = sumdif / clusternum;%�в�
   threshold = 0.01;%��ֵ
   if res < threshold
        break;
   end
   last_center_points = center_points;
   %��ʾ�м����
   tempborder = boundarymask(cluster_array);
   ima_show = ima;%չʾ���м����ͼ��
   %�۵����ģ���ɫ���
   for i = 1:clusternum
        x = center_points(i,1);
        y = center_points(i,2);
        roundpoints = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y;x,y+1;x+1,y-1;x+1,y;x+1,y+1];%3*3��Χ�ڵĵ�
        for k = 1:9
            c_point = correct_points(roundpoints(k,:),height,width);
            ima_show(c_point(1),c_point(2),:) = [255,0,0];
        end
   end
   %�߽磬��ɫ���
   for i = 1:height
       for j = 1:width
           if tempborder(i,j) == 1
                ima_show(i,j,:) = [83,117,213];
           end
       end
   end
   figure,imshow(ima_show);
end
%%
%����Ż�
squ = strel('square',8);
for i =1:clusternum
    mask = zeros(height,width);
    index1 = find(cluster_array == i);
    mask(index1) = 1;%ÿһ���mask
    mask2 = imclose(mask,squ);%�ȱ�
    mask3 = imopen(mask2,squ);%�ٿ�
    index2 = find(mask3);
    cluster_array(index2) = i;
end
%%
%���
cluster = cluster_array;
border = boundarymask(cluster_array);
end