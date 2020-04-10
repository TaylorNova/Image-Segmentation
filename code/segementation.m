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
axes(handles.axes1) %��TagֵΪaxes1����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
axes(handles.axes2) %��TagֵΪaxes2����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
axes(handles.axes3) %��TagֵΪaxes3����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
axes(handles.axes4) %��TagֵΪaxes4����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
axes(handles.axes5) %��TagֵΪaxes5����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
axes(handles.axes7) %��TagֵΪaxes7����������Ϊ��ǰ
imshow(initial_ima); %��ʾͼƬ
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
temppicture = get(handles.chooselist,'value');%ѡ���ͼƬ
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
%mask��ʼ��
maskfro = zeros(height,width);
maskback = zeros(height,width);
%��ʾͼƬ
axes(handles.axes1) %��TagֵΪaxes1����������Ϊ��ǰ
imshow(ima);
clusternum = 1000;%��������
if temppicture == 3
    clusternum = clusternum + 15;
end
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
   axes(handles.axes7) %��TagֵΪaxes7����������Ϊ��ǰ
   imshow(ima_show);
end
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
%���
border = boundarymask(cluster_array);
%��ʾ�����ؽ��
outputImage = zeros(size(ima),'like',ima);
for labelVal = 1:clusternum
    redIdx = find(cluster_array == labelVal);
    greenIdx = redIdx+height*width;
    blueIdx = redIdx+2*height*width;
    outputImage(redIdx) = mean(ima(redIdx));
    outputImage(greenIdx) = mean(ima(greenIdx));
    outputImage(blueIdx) = mean(ima(blueIdx));
end
axes(handles.axes3) %��TagֵΪaxes3����������Ϊ��ǰ
imshow(outputImage,'InitialMagnification',67);
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in choosefro.
function choosefro_Callback(hObject, eventdata, handles)
% hObject    handle to choosefro (see GCBO)
global height width imastr;
global maskfro;
ima = imread(imastr);
axes(handles.axes1) %��TagֵΪaxes1����������Ϊ��ǰ
imshow(ima); %��ʾͼƬ
tempmask = zeros(height,width);
%ѡ���������
h = imrect;
pos = wait(h);
pos = round(pos);
tempmask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3)) = 1;
maskfro = maskfro | tempmask;
axes(handles.axes4) %��TagֵΪaxes4����������Ϊ��ǰ
imshow(maskfro); %��ʾͼƬ
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chooseback.
function chooseback_Callback(hObject, eventdata, handles)
% hObject    handle to chooseback (see GCBO)
global height width imastr;
global maskback;
ima = imread(imastr);
axes(handles.axes1) %��TagֵΪaxes1����������Ϊ��ǰ
imshow(ima); %��ʾͼƬ
tempmask = zeros(height,width);
%ѡ���������
h = imrect;
pos = wait(h);
pos = round(pos);
tempmask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3)) = 1;
maskback = maskback | tempmask;
axes(handles.axes5) %��TagֵΪaxes5����������Ϊ��ǰ
imshow(maskback); %��ʾͼƬ
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
global height width imastr;
global maskback maskfro;
global cluster_array;
ima = imread(imastr);%ԭͼ
ima_show = ima;%չʾͼ
BW = lazysnapping(ima,cluster_array,maskfro,maskback);%���˿�ͼ
for i = 1:height
    for j = 1:width
        if BW(i,j) == 0
            ima_show(i,j,:) = [0,0,0];
        end
    end
end
axes(handles.axes2) %��TagֵΪaxes2����������Ϊ��ǰ
imshow(ima_show); %��ʾͼƬ
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
