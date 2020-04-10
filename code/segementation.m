function varargout = segementation(varargin)
% SEGEMENTATION MATLAB code for segementation.fig
%      SEGEMENTATION, by itself, creates a new SEGEMENTATION or raises the existing
%      singleton*.
%
%      H = SEGEMENTATION returns the handle to a new SEGEMENTATION or the handle to
%      the existing singleton*.
%
%      SEGEMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGEMENTATION.M with the given input arguments.
%
%      SEGEMENTATION('Property','Value',...) creates a new SEGEMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segementation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segementation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  chooseimage "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segementation

% Last Modified by GUIDE v2.5 11-Dec-2019 22:37:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segementation_OpeningFcn, ...
                   'gui_OutputFcn',  @segementation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before segementation is made visible.
function segementation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segementation (see VARARGIN)

% chooseimage default command line output for segementation
handles.output = hObject;
global height width imastr;
global maskfro maskback;
global cluster_array;
initial_ima = imread('0.jpg');
axes(handles.axes1) %将Tag值为axes1的坐标轴置为当前
imshow(initial_ima); %显示图片
axes(handles.axes2) %将Tag值为axes2的坐标轴置为当前
imshow(initial_ima); %显示图片
axes(handles.axes3) %将Tag值为axes3的坐标轴置为当前
imshow(initial_ima); %显示图片
axes(handles.axes4) %将Tag值为axes4的坐标轴置为当前
imshow(initial_ima); %显示图片
axes(handles.axes5) %将Tag值为axes5的坐标轴置为当前
imshow(initial_ima); %显示图片
axes(handles.axes7) %将Tag值为axes7的坐标轴置为当前
imshow(initial_ima); %显示图片
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segementation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = segementation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in chooselist.
function chooselist_Callback(hObject, eventdata, handles)
% hObject    handle to chooselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chooselist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chooselist


% --- Executes during object creation, after setting all properties.
function chooselist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooselist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chooseimage.
function chooseimage_Callback(hObject, eventdata, handles)
% hObject    handle to chooseimage (see GCBO)
temppicture = get(handles.chooselist,'value');%选择的图片
if temppicture == 1
    str = '.\data\1.jpg';
elseif temppicture == 2
    str = '.\data\2.jpg';
elseif temppicture == 3
    str = '.\data\3.jpg';
end
global height width imastr;
global cluster_array;
global maskfro maskback;
imastr = str;
ima = imread(str);
[height,width,channel] = size(ima);
%mask初始化
maskfro = zeros(height,width);
maskback = zeros(height,width);
%显示图片
axes(handles.axes1) %将Tag值为axes1的坐标轴置为当前
imshow(ima);
clusternum = 1000;%聚类数量
if temppicture == 3
    clusternum = clusternum + 15;
end
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
   axes(handles.axes7) %将Tag值为axes7的坐标轴置为当前
   imshow(ima_show);
end
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
%结果
border = boundarymask(cluster_array);
%显示超像素结果
outputImage = zeros(size(ima),'like',ima);
for labelVal = 1:clusternum
    redIdx = find(cluster_array == labelVal);
    greenIdx = redIdx+height*width;
    blueIdx = redIdx+2*height*width;
    outputImage(redIdx) = mean(ima(redIdx));
    outputImage(greenIdx) = mean(ima(greenIdx));
    outputImage(blueIdx) = mean(ima(blueIdx));
end
axes(handles.axes3) %将Tag值为axes3的坐标轴置为当前
imshow(outputImage,'InitialMagnification',67);
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in choosefro.
function choosefro_Callback(hObject, eventdata, handles)
% hObject    handle to choosefro (see GCBO)
global height width imastr;
global maskfro;
ima = imread(imastr);
axes(handles.axes1) %将Tag值为axes1的坐标轴置为当前
imshow(ima); %显示图片
tempmask = zeros(height,width);
%选择矩形区域
h = imrect;
pos = wait(h);
pos = round(pos);
tempmask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3)) = 1;
maskfro = maskfro | tempmask;
axes(handles.axes4) %将Tag值为axes4的坐标轴置为当前
imshow(maskfro); %显示图片
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chooseback.
function chooseback_Callback(hObject, eventdata, handles)
% hObject    handle to chooseback (see GCBO)
global height width imastr;
global maskback;
ima = imread(imastr);
axes(handles.axes1) %将Tag值为axes1的坐标轴置为当前
imshow(ima); %显示图片
tempmask = zeros(height,width);
%选择矩形区域
h = imrect;
pos = wait(h);
pos = round(pos);
tempmask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3)) = 1;
maskback = maskback | tempmask;
axes(handles.axes5) %将Tag值为axes5的坐标轴置为当前
imshow(maskback); %显示图片
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
global height width imastr;
global maskback maskfro;
global cluster_array;
ima = imread(imastr);%原图
ima_show = ima;%展示图
BW = lazysnapping(ima,cluster_array,maskfro,maskback);%懒人抠图
for i = 1:height
    for j = 1:width
        if BW(i,j) == 0
            ima_show(i,j,:) = [0,0,0];
        end
    end
end
axes(handles.axes2) %将Tag值为axes2的坐标轴置为当前
imshow(ima_show); %显示图片
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
