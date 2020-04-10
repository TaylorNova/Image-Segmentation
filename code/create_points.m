function points = create_points(height,width,clusternum)
%����ͼ���Ͼ��ȷֲ��ĵ�ĺ���
    S = (height*width/clusternum)^0.5;%Sֵ�ļ���
    m_step = round(height/S);%����
    n_step = round(width/S);
    points = [];%��ʼ����
    startpoint = [S/2,S/2];%��ʼ��
    for i = 1:m_step
        for j = 1:n_step
            temppoints = [round(startpoint(1) + (i-1)*S),round(startpoint(2) + (j-1)*S)];
            points = [points;temppoints];
        end
    end
    shape = size(points);
    len = shape(1);%�Ѿ������ĵ���
    if len < clusternum
        randlist = randperm(len);%�����
        %����ʣ�µĵ�
        for i = 1:clusternum-len
            temppoint = [round(points(randlist(i),1)-S/3),round(points(randlist(i),2)-S/3)];
            points = [points;temppoint];
        end
    end
end