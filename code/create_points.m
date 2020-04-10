function points = create_points(height,width,clusternum)
%生成图像上均匀分布的点的函数
    S = (height*width/clusternum)^0.5;%S值的计算
    m_step = round(height/S);%步长
    n_step = round(width/S);
    points = [];%初始化点
    startpoint = [S/2,S/2];%起始点
    for i = 1:m_step
        for j = 1:n_step
            temppoints = [round(startpoint(1) + (i-1)*S),round(startpoint(2) + (j-1)*S)];
            points = [points;temppoints];
        end
    end
    shape = size(points);
    len = shape(1);%已经产生的点数
    if len < clusternum
        randlist = randperm(len);%随机数
        %生成剩下的点
        for i = 1:clusternum-len
            temppoint = [round(points(randlist(i),1)-S/3),round(points(randlist(i),2)-S/3)];
            points = [points;temppoint];
        end
    end
end