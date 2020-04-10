function [cluster,border] = mySLIC(ima,clusternum)
%%
%预处理
I = ima;%原图
I = im2double(I);
I_lab = rgb2lab(I);%LAB图像
[height,width,channel] = size(ima);
mvalue = 20;%m值
S = round((height*width/clusternum)^0.5);%S值的计算
%Gabor滤波
wavelength = 2.^(3:5) * 3;
orientation = 0:90:90;
g = gabor(wavelength,orientation); %2 * 3 = 6个
I = rgb2gray(im2single(ima));
gabormag = imgaborfilt(I,g);
for i = 1:6
    maxvalue = max(max(gabormag(:,:,i)));
    minvalue = min(min(gabormag(:,:,i)));
    len = maxvalue - minvalue;
    gabormag(:,:,i) = (gabormag(:,:,i) - minvalue)/len;%归一化
end
%%
%聚类中心点初始化
center_points = create_points(height,width,clusternum);
%初始点梯度纠偏
I_grad = gradient(I_lab);%梯度
I_grad_abs = abs(I_grad);%梯度幅度
for i = 1:clusternum
    x = center_points(i,1);
    y = center_points(i,2);
    roundpoints = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y;x,y+1;x+1,y-1;x+1,y;x+1,y+1];%3*3范围内的点
    mingrad = 10000000;%范围内最小梯度
    for k = 1:9
        c_point = correct_points(roundpoints(k,:),height,width);
        if I_grad_abs(c_point(1),c_point(2)) < mingrad
            mingrad = I_grad_abs(c_point(1),c_point(2));
        end
    end
    for k = 1:9 %最小梯度对应点
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
%循环聚类
last_center_points = center_points;%上一次的聚类中心点
cluster_array = -1*ones(height,width);%聚类标签初始化为-1
distance_array = inf*ones(height,width);%距离初始化为无穷
while true
    %计算距离，更新聚类标签
    for i = 1:clusternum
        %x，y的搜索范围
        x_range = [max(center_points(i,1) - S,1),min(center_points(i,1) + S,height)];
        y_range = [max(center_points(i,2) - S,1),min(center_points(i,2) + S,width)];
        for m = x_range(1):x_range(2)
            for n = y_range(1):y_range(2)
                centerlab = I_lab(center_points(i,1),center_points(i,2),:);%中心点的lab值
                templab = I_lab(m,n,:);%当前点的lab值
                lab_dis = (centerlab(1) - templab(1))^2 + (centerlab(2) - templab(2))^2 + (centerlab(3) - templab(3))^2;%lab距离的平方
                space_dis = (center_points(i,1) - m)^2 + (center_points(i,2) - n)^2;%空间距离的平方
                sumdis = 0;
                for v = 1:6
                    sumdis = sumdis + (gabormag(m,n,v) - gabormag(center_points(i,1),center_points(i,2),v))^2;
                end
                gabor_dis = sumdis/6;%gabor滤波距离的平方
                total_dis = (lab_dis/(mvalue^2) + space_dis/(S^2) + gabor_dis)^0.5;%总距离
                if total_dis < distance_array(m,n)
                    distance_array(m,n) = total_dis;
                    cluster_array(m,n) = i;
                end
            end
        end
    end
   %更新距离中心点
   count_sum = zeros(clusternum,2); %各类点x，y坐标总和
   countnum = zeros(clusternum,1); %各类点总数
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
   %均值法更新聚类中心点
   for i = 1:clusternum
        center_points(i,1) = round(count_sum(i,1)/countnum(i));
        center_points(i,2) = round(count_sum(i,2)/countnum(i));
   end
   %计算残差，判断是否结束循环
   dif = center_points - last_center_points;
   sumdif = 0;
   for i = 1:clusternum
       sumdif = sumdif + (dif(i,1)^2 + dif(i,2)^2)^0.5;
   end
   res = sumdif / clusternum;%残差
   threshold = 0.01;%阈值
   if res < threshold
        break;
   end
   last_center_points = center_points;
   %显示中间过程
   tempborder = boundarymask(cluster_array);
   ima_show = ima;%展示的中间过程图像
   %聚点中心，红色标记
   for i = 1:clusternum
        x = center_points(i,1);
        y = center_points(i,2);
        roundpoints = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y;x,y+1;x+1,y-1;x+1,y;x+1,y+1];%3*3范围内的点
        for k = 1:9
            c_point = correct_points(roundpoints(k,:),height,width);
            ima_show(c_point(1),c_point(2),:) = [255,0,0];
        end
   end
   %边界，蓝色标记
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
%结果优化
squ = strel('square',8);
for i =1:clusternum
    mask = zeros(height,width);
    index1 = find(cluster_array == i);
    mask(index1) = 1;%每一类的mask
    mask2 = imclose(mask,squ);%先闭
    mask3 = imopen(mask2,squ);%再开
    index2 = find(mask3);
    cluster_array(index2) = i;
end
%%
%结果
cluster = cluster_array;
border = boundarymask(cluster_array);
end