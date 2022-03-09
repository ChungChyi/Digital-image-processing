function varargout = DGP(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DGP_OpeningFcn, ...
                   'gui_OutputFcn',  @DGP_OutputFcn, ...
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
function DGP_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);                       
axes(handles.pic1);                                        %初始化原始图像
imshow(uint8(255));
axes(handles.pic2);                                    %初始化噪声图像
imshow(255*ones(255,255));
axes(handles.pic3);                                       %初始化处理图像
imshow(255*ones(255,255));
set(handles.NoiseChoose,'Enable','off');                          %图片未选择时，滤波方法，添加噪声选项，对比方法不可选
set(handles.CoChange,'Enable','off');
set(handles.HisSta,'Enable','off');
set(handles.FreFiltering,'Enable','off');
set(handles.SpaFiltering,'Enable','off');
set(handles.ColorEn,'Enable','off');
global Size ProNimg Proimg img val;
Size=3;ProNimg=[255*ones(255,255)];            %全局变量的初始化，Size为模板大小，ProNimg为添加噪声后图像
img=[255*ones(255,255)];                                 %img为选择图片，Proimg为滤波后图像
Proimg=[255*ones(255,255)];
val=0;   

function varargout = DGP_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
% --------------------------------------------------------------------
%%文件菜单
function file_Callback(hObject, eventdata, handles)
%文件-打开函数
function Open_Callback(hObject, eventdata, handles)
global img;
[name,path]=uigetfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                '选择图片');                    %读取图像
img=imread([path,name]);
axes(handles.pic1);
imshow(img);                                                      %显示图像并开启不可选菜单
set(handles.NoiseChoose,'Enable','on');                          %图片未选择时，滤波方法，添加噪声选项，对比方法不可选
set(handles.CoChange,'Enable','on');
set(handles.HisSta,'Enable','on');
set(handles.FreFiltering,'Enable','on');
set(handles.SpaFiltering,'Enable','on');
set(handles.ColorEn,'Enable','on');
%文件-退出函数
function Quit_Callback(hObject, eventdata, handles)
close(gcf);

% --------------------------------------------------------------------
%%噪声选择菜单
function NoiseChoose_Callback(hObject, eventdata, handles)   
%%噪声选择-高斯噪声 菜单
function GN_Callback(hObject, eventdata, handles)
   %%默认方差0.01 函数
   function GMFC_Callback(hObject, eventdata, handles)
global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'gaussian',0,0.01));
axes(handles.pic2);
imshow(ProNimg);   
   %%自选参数 函数
   function GAC_Callback(hObject, eventdata, handles)
     global img ProNimg val;
val=1;
str=['请输入高斯噪声方差'];
title=['高斯噪声自选参数'];
input=inputdlg(str,title,1,{'0.01'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'gaussian',0,input));
axes(handles.pic2);
imshow(ProNimg);
%%噪声选择-随机噪声 菜单
function RN_Callback(hObject, eventdata, handles)
    %%默认方差0.01 函数
    function RMFC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'speckle',0.01));
axes(handles.pic2);
imshow(ProNimg);    
    %%自选参数  函数
    function RAC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
str=['请输入随机噪声方差'];
title=['随机噪声噪声自选参数'];
input=inputdlg(str,title,1,{'0.04'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'speckle',input));
axes(handles.pic2);
imshow(ProNimg);
%%噪声选择-椒盐噪声 菜单
function SN_Callback(hObject, eventdata, handles)
    %%默认强度0.05 函数
    function MS_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'salt & pepper',0.05));
axes(handles.pic2);
imshow(ProNimg);
    %%自选参数 函数
    function SAC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
str=['请输入椒盐噪声强度'];
title=['椒盐噪声自选参数'];
input=inputdlg(str,title,1,{'0.05'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'salt & pepper',input));
axes(handles.pic2);
imshow(ProNimg);    
%%噪声选择-泊松噪声 函数
    function PN_Callback(hObject, eventdata, handles)
global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'poisson'));
axes(handles.pic2);
imshow(ProNimg);    



% --------------------------------------------------------------------
%%坐标变换菜单
function CoChange_Callback(hObject, eventdata, handles)
  %%坐标变换-灰度变换 函数
  function Graychange_Callback(hObject, eventdata, handles)
  global img   ProNimg val
  if size(img,3)==3
  ProNimg=rgb2gray(img);
  val=1;
  axes(handles.pic2);
  imshow(ProNimg);
  xlabel('gray image');
  else
      msgbox('is gray image','fail');
  end
  %%
  function LCF_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
  axes(handles.pic2);
prompt={'Input magnification:'};
defans={'2'};
p=inputdlg(prompt,'magnification',1,defans);
p1=str2num(p{1});
ProNimg=imresize(img,p1,'nearest');  
imshow(ProNimg);
xlabel('nearest neighbor');    
  %%
  function DLF_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
      axes(handles.pic2);
prompt={'Input magnification:'};
defans={'2'};
p=inputdlg(prompt,'magnification',1,defans);
p1=str2num(p{1});
ProNimg=imresize(img,p1,'bilinear');           
imshow(ProNimg);
xlabel('Bilinear Interpolation');
  %%
  function DTF_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
  axes(handles.pic2);
  prompt={'Input magnification:'};
defans={'2'};
p=inputdlg(prompt,'magnification',1,defans);
p1=str2num(p{1});
ProNimg=imresize(img,p1,'bicubic');           
imshow(ProNimg);
xlabel('Bicubic interpolation');
  %%
  function UPD_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
  axes(handles.pic2);
if size(img,3)==3 
for k=1:3 
y(:,:,k)=flipud(img(:,:,k));
end
imshow(y);
xlabel('Upside down');
else
    ProNimg=flipud(img); 
    imshow(ProNimg);
    xlabel('Upside down');
end    
  %%
  function LRV_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
      axes(handles.pic2);
  if size(img,3)==3
for k=1:3 
y(:,:,k)=fliplr(img(:,:,k));
end
imshow(y);
xlabel('Turn around');
  else
    ProNimg=fliplr(img);
    imshow(ProNimg);
    xlabel('Turn around');
end    
  %%
  function AAO_Callback(hObject, eventdata, handles)
      global img val ProNimg
      val=1;
  axes(handles.pic2);
prompt={'input parameter 1:'};
defans={'30'};
p=inputdlg(prompt,'input parameter',1,defans);
p1=str2num(p{1});
ProNimg=imrotate(img,p1);          
imshow(ProNimg);
xlabel(p1);

% --------------------------------------------------------------------
%%直方图统计菜单
function HisSta_Callback(hObject, eventdata, handles)
    %%直方图统-R直方图 函数
    function RHis_Callback(hObject, eventdata, handles)
        global img
set(handles.pic3,'HandleVisibility','ON');
axes(handles.pic3);
R=img(:,:,1);
    x=imhist(R);                       
    x1=x(1:10:256);
   horz=1:10:256;
   bar(horz,x1);
   xlabel('R histogram');
   set(handles.pic3,'xtick',0:50:255);
    %%直方图统-G直方图 函数
    function GHis_Callback(hObject, eventdata, handles)
        global img
set(handles.pic3,'HandleVisibility','ON');
axes(handles.pic3);
G=img(:,:,2);
    x=imhist(G);                       
    x1=x(1:10:256);
   horz=1:10:256;
   bar(horz,x1);
   xlabel('G histogram');
   set(handles.pic3,'xtick',0:50:255);
    %%直方图统-B直方图 函数
    function BHis_Callback(hObject, eventdata, handles)
        global img
set(handles.pic3,'HandleVisibility','ON');
axes(handles.pic3);
B=img(:,:,3);
    x=imhist(B);                       
    x1=x(1:10:256);
   horz=1:10:256;
   bar(horz,x1);
   xlabel('B histogram');
   set(handles.pic3,'xtick',0:50:255);
    %%直方图统-直方图 函数
    function His_Callback(hObject, eventdata, handles)
        global img ProNimg val
        if val==0
            ProNimg=img;
        end
        axes(handles.pic3);
        if size(ProNimg,3)==1
        H=histeq(ProNimg);
        else H(:,:,1)=histeq(ProNimg(:,:,1));
            H(:,:,2)=histeq(ProNimg(:,:,2));
            H(:,:,3)=histeq(ProNimg(:,:,3));
        end
        imshow(H);

        
% --------------------------------------------------------------------
%%频率域滤波增强
function FreFiltering_Callback(hObject, eventdata, handles)
    %%频域滤波-频谱图 函数
    function FSP_Callback(hObject, eventdata, handles)
        global img
        axes(handles.pic3);
if size(img,3)==3
    m=fft2(img(:,:,1));
    y=fftshift(m);
    imshow(log(abs(y)),[]);
    xlabel('Spectrum');
else
    m=fft2(img);
    y=fftshift(m);
    imshow(log(abs(y)),[]);
     xlabel('Spectrum');
end
    %%频域滤波-巴特沃斯低通滤波 函数
    function BLF_Callback(hObject, eventdata, handles)
        global img ProNimg val
    axes(handles.pic3);
    if val==0
        ProNimg=img;
    end
        x=rgb2gray(ProNimg);
        f=double(x);
        g=fft2(f);
        g=fftshift(g);
        [M,N]=size(g);
        nn=2;
        d0=50;
        m=fix(M/2); n=fix(N/2);
        for i=1:M
            for j=1:N  
                d=sqrt((i-m)^2+(j-n)^2);
                 h=1/(1+0.414*(d/d0)^(2*nn));
                result(i,j)=h*g(i,j);
            end 
        end
    result=ifftshift(result);
    y2=ifft2(result);
    Proimg=uint8(real(y2));
    imshow(Proimg);
    xlabel('Butterworth low-pass filter');    
    %%频域滤波-巴特沃斯高通滤波 函数
    function BHF_Callback(hObject, eventdata, handles)
        global img ProNimg val
    axes(handles.pic3);
    if val==0
        ProNimg=img;
    end
        x=rgb2gray(ProNimg);
        f=double(x);
        g=fft2(f);
        g=fftshift(g);
        [M,N]=size(g);
        nn=2;
        d0=50;
        m=fix(M/2); n=fix(N/2);
        for i=1:M
            for j=1:N  
                d=sqrt((i-m)^2+(j-n)^2);
                 if d<=d0
                h=0;
            else h=1;                
            end
                result(i,j)=h*g(i,j);
            end 
        end
    result=ifftshift(result);
    y2=ifft2(result);
    Proimg=uint8(real(y2));
    imshow(Proimg);
    xlabel('Butterworth high-pass filter');
    %%频域滤波-离散傅里叶变换 函数
    function DFF_Callback(hObject, eventdata, handles)
        global img 
    axes(handles.pic3);
    i=im2double(img);
    f1=fft2(i);
    fc1=fftshift(f1);
    i=log(1+abs(fc1));
    imshow(i);
    xlabel('The Fourier transform');
    %%频域滤波-离散余弦变换 函数
    function DCC_Callback(hObject, eventdata, handles)
        global img 
    axes(handles.pic3);
    if size(img,3)==1
        msgbox('It is a gray image!','Error');
    elseif size(img,3)==3
        i=rgb2gray(img);
        d=dct2(i);
        i=log(abs(d));
        imshow(i);
        xlabel('Discrete cosine transform');
        else 
    msgbox('It is not suitable source image or source image have been dealing with','Error');
end
        
 
% --------------------------------------------------------------------
%%空间域滤波增强 菜单
function SpaFiltering_Callback(hObject, eventdata, handles)
%%空间域滤波-平滑滤波方法 菜单
function PF_Callback(hObject, eventdata, handles)
    %%空间域滤波-平滑滤波方法-均值滤波 函数
    function UF_Callback(hObject, eventdata, handles)
        global ProNimg img  val;
 axes(handles.pic3)
    if val==0;
        ProNimg=img;
    end
    prompt={'Input ModeSize:'};
defans={'5'};
p=inputdlg(prompt,'ModeSize',1,defans);
Size=str2num(p{1});
filtermode=fspecial('average',[Size,Size]);
if size(ProNimg,3)==3
       Proimg(:,:,1)=filter2(filtermode,ProNimg(:,:,1),'same');
       Proimg(:,:,2)=filter2(filtermode,ProNimg(:,:,2),'same');
       Proimg(:,:,3)=filter2(filtermode,ProNimg(:,:,3),'same');
        imshow(Proimg);
    else Proimg=filter2(filtermode,ProNimg,'same');
        imshow(Proimg);
end
Proimg=uint8(Proimg);
imshow(Proimg);
    %%空间域滤波-平滑滤波方法-中值滤波 函数
    function MF_Callback(hObject, eventdata, handles)
    global img ProNimg val
    axes(handles.pic3)
    if val==0;
        ProNimg=img;
    end
    prompt={'Input ModeSize:'};
defans={'5'};
p=inputdlg(prompt,'ModeSize',1,defans);
Size=str2num(p{1});
    if size(ProNimg,3)==3
        r=medfilt2(ProNimg(:,:,1),[Size,Size]);
        g=medfilt2(ProNimg(:,:,2),[Size,Size]);
        b=medfilt2(ProNimg(:,:,3),[Size,Size]);
       Proimg(:,:,1)=r;
       Proimg(:,:,2)=g;
       Proimg(:,:,3)=b;
        imshow(Proimg);
    else Proimg=medfilt2(ProNimg,[Size,Size]);
        imshow(Proimg);
    end
    %%空间域滤波-平滑滤波方法-自适应滤波 函数
    function ZF_Callback(hObject, eventdata, handles)
    global img ProNimg val
    axes(handles.pic3)
    if val==0;
        ProNimg=img;
    end
    prompt={'Input ModeSize:'};
defans={'3'};
p=inputdlg(prompt,'ModeSize',1,defans);
Size=str2num(p{1});
    if size(ProNimg,3)==3
        Proimg(:,:,1)=wiener2(ProNimg(:,:,1),[Size,Size]);
        Proimg(:,:,2)=wiener2(ProNimg(:,:,2),[Size,Size]);
        Proimg(:,:,3)=wiener2(ProNimg(:,:,3),[Size,Size]);
        imshow(Proimg);
    else Proimg=wiener2(ProNimg,[Size,Size]);
        imshow(Proimg);
    end    
%%空间域滤波-锐化滤波 菜单
function RF_Callback(hObject, eventdata, handles)
    %%空间域滤波-锐化滤波-梯度滤波 函数
    function TF_Callback(hObject, eventdata, handles)
  global img ProNimg Proimg val;
if val==0
    ProNimg=img;
end

if size(ProNimg,3)==3
[a b]=size(ProNimg(:,:,1));
        im=ProNimg(:,:,1);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,1)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
        [a b]=size(ProNimg(:,:,2));
        im=ProNimg(:,:,2);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,2)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
        [a b]=size(ProNimg(:,:,3));
        im=ProNimg(:,:,3);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,3)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
else
    [a b]=size(ProNimg);
        im=ProNimg;                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
end
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg); 
    %%空间域滤波-锐化滤波-梯度锐化 函数
    function TR_Callback(hObject, eventdata, handles)
        global img ProNimg Proimg val;
if val==0
    ProNimg=img;
end
im1=ProNimg;
if size(ProNimg,3)==3
[a b]=size(ProNimg(:,:,1));
        im=ProNimg(:,:,1);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,1)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
        [a b]=size(ProNimg(:,:,2));
        im=ProNimg(:,:,2);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,2)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
        [a b]=size(ProNimg(:,:,3));
        im=ProNimg(:,:,3);                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j,3)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
else
    [a b]=size(ProNimg);
        im=ProNimg;                                              
        for i=1:a                                                  
            for j=1:b
                if (i+1)==(a+1)
                    im(i+1,j+1)=0;im(i+1,j)=0;
                else if (j+1)==(b+1)
                        im(i+1,j+1)=0;im(i,j+1)=0;
                    end
                end
                Proimg(i,j)=abs(im(i,j)-im(i+1,j+1))+abs(im(i+1,j)-im(i,j+1));
            end
        end
end
Proimg=double(Proimg)+double(im1);
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg);
    %%空间域滤波-锐化滤波-拉普拉斯算子滤波 函数
    function LF_Callback(hObject, eventdata, handles)
    global img ProNimg val;
if val==0
    ProNimg=img;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
end
if size(ProNimg,3)==3
filtermodel=[0 -1 0;-1 4 -1;0 -1 0];
Proimg(:,:,1)=filter2(filtermodel,ProNimg(:,:,1),'same');
Proimg(:,:,2)=filter2(filtermodel,ProNimg(:,:,2),'same');
Proimg(:,:,3)=filter2(filtermodel,ProNimg(:,:,3),'same');
else
    Proimg=filter2(filtermodel,ProNimg,'same');
end
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg);    
    %%空间域滤波-锐化滤波-拉普拉斯算子锐化 函数 
    function LR_Callback(hObject, eventdata, handles)
    global img ProNimg val;
if val==0
    ProNimg=img;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
end
im=ProNimg;
filtermodel=[0 -1 0;-1 4 -1;0 -1 0];
if size(ProNimg,3)==3
Proimg(:,:,1)=filter2(filtermodel,ProNimg(:,:,1),'same');
Proimg(:,:,2)=filter2(filtermodel,ProNimg(:,:,2),'same');
Proimg(:,:,3)=filter2(filtermodel,ProNimg(:,:,3),'same');
else
    Proimg=filter2(filtermodel,ProNimg,'same');
end
Proimg=Proimg+double(im);
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg)
    %%空间域滤波-锐化滤波-定向滤波 菜单
    function DF_Callback(hObject, eventdata, handles)
        %%空间域滤波-锐化滤波-定向滤波-垂直定向滤波 函数
        function VDF_Callback(hObject, eventdata, handles)
         global img ProNimg  val;
if val==0
    ProNimg=img;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
end
im=ProNimg;
filtermodel=[-1 -1 -1;2 2 2;-1 -1 -1];
if size(ProNimg,3)==3
Proimg(:,:,1)=filter2(filtermodel,ProNimg(:,:,1),'same');
Proimg(:,:,2)=filter2(filtermodel,ProNimg(:,:,2),'same');
Proimg(:,:,3)=filter2(filtermodel,ProNimg(:,:,3),'same');
else
    Proimg=filter2(filtermodel,ProNimg,'same');
end
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg);   
        %%空间域滤波-锐化滤波-定向滤波-水平定向滤波 函数
        function LDF_Callback(hObject, eventdata, handles)
            global img ProNimg  val;
if val==0
    ProNimg=img;
end
filtermodel=[-1 2 -1;-1 2 -1;-1 2 -1];
if size(ProNimg,3)==3
Proimg(:,:,1)=filter2(filtermodel,ProNimg(:,:,1),'same');
Proimg(:,:,2)=filter2(filtermodel,ProNimg(:,:,2),'same');
Proimg(:,:,3)=filter2(filtermodel,ProNimg(:,:,3),'same');
else
    Proimg=filter2(filtermodel,ProNimg,'same');
end
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg);
        %%空间域滤波-锐化滤波-定向滤波-对角定向滤波 函数    
        function DDF_Callback(hObject, eventdata, handles)
            global img ProNimg  val;
if val==0
    ProNimg=img;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
end
filtermodel=[-1 -1 2;-1 2 -1;2 -1 -1];
if size(ProNimg,3)==3
Proimg(:,:,1)=filter2(filtermodel,ProNimg(:,:,1),'same');
Proimg(:,:,2)=filter2(filtermodel,ProNimg(:,:,2),'same');
Proimg(:,:,3)=filter2(filtermodel,ProNimg(:,:,3),'same');
else
    Proimg=filter2(filtermodel,ProNimg,'same');
end
Proimg=uint8(Proimg);
axes(handles.pic3);
imshow(Proimg);
        
  


% --------------------------------------------------------------------
%%彩色增强菜单
function ColorEn_Callback(hObject, eventdata, handles)
    %%彩色增强-真彩色增强 菜单
    function TCH_Callback(hObject, eventdata, handles)
        %%彩色增强-真彩色增强-红色分量 函数
        function TRed_Callback(hObject, eventdata, handles)
        global img   
axes(handles.pic3);
imshow(img(:,:,1));
xlabel('red component');
        %%彩色增强-真彩色增强-蓝色分量 函数    
        function TBlue_Callback(hObject, eventdata, handles)
         global img       
axes(handles.pic3);
imshow(img(:,:,2));
xlabel('blue component');
        %%彩色增强-真彩色增强-绿色分量 函数
        function TGreen_Callback(hObject, eventdata, handles)
                   global img       
axes(handles.pic3);
imshow(img(:,:,3));
xlabel('green component');
%%彩色增强-伪彩色增强 菜单
function FCH_Callback(hObject, eventdata, handles)
    %%彩色增强-伪彩色增-密度分割 函数
    function DS_Callback(hObject, eventdata, handles)
global img ProNimg Proimg
ProNimg=rgb2gray(img);
[m,n]=size(ProNimg);
r=zeros(m,n);
g=zeros(m,n);
b=zeros(m,n);
e=max(max(ProNimg));
f=min(min(ProNimg));
t=(e-f)/3;
for i=1:m
    for j=1:n
        if ProNimg(i,j)<(f+t);
            r(i,j)=255;
        else if ProNimg(i,j)>(e-t)
                g(i,j)=255;
            else
                b(i,j)=255;
            end
        end
    end
end
Proimg=uint8(zeros(m,n,3));
Proimg(:,:,1)=uint8(r);
Proimg(:,:,2)=uint8(g);
Proimg(:,:,3)=uint8(b);
axes(handles.pic3);
imshow(Proimg);
xlabel('Density slicing');
    %%彩色增强-伪彩色增-灰度级 函数
    function HDF_Callback(hObject, eventdata, handles)
    global img Proimg
    w=rgb2gray(img);
    [m,n]=size(w);
    r=w;
    g=w;
    b=w;
for i=1:m
    for j=1:n
        if w(i,j)<=127
            r(i,j)=0;
        else if w(i,j)>=191
                r(i,j)=255;
            else
                r(i,j)=(w(i,j)-127)*4;
            end
        end
    end
end
for i=1:m
    for j=1:n
        if w(i,j)<=63
            g(i,j)=4*w(i,j);
        else if w(i,j)>=192
                g(i,j)=(255-w(i,j))*4;
            else
                g(i,j)=255;
            end
        end
    end
end
for i=1:m
    for j=1:n
        if w(i,j)<=63
            b(i,j)=255;
        else if w(i,j)>=127
                b(i,j)=0;
            else
                b(i,j)=(127-w(i,j))*4;
            end
        end
    end
end
Proimg=uint8(zeros(m,n,3));
Proimg(:,:,1)=uint8(r);
Proimg(:,:,2)=uint8(g);
Proimg(:,:,3)=uint8(b);
axes(handles.pic3);
imshow(Proimg);
xlabel('Density slicing');
    
    
% 选择图片按钮
function choose_Callback(hObject, eventdata, handles)
    global img;
[name,path]=uigetfile({'*.jpg';'*.bmp'},'选择图片');                    %读取图像
img=imread([path,name]);
axes(handles.pic1);
imshow(img);                                                      %显示图像并开启不可选菜单
set(handles.NoiseChoose,'Enable','on');                          %图片未选择时，滤波方法，添加噪声选项，对比方法不可选
set(handles.CoChange,'Enable','on');
set(handles.HisSta,'Enable','on');
set(handles.FreFiltering,'Enable','on');
set(handles.SpaFiltering,'Enable','on');
set(handles.ColorEn,'Enable','on');
%[图片变换]清除图片按钮
function clear1_Callback(hObject, eventdata, handles)
    global val
    val=0;
axes(handles.pic2);
imshow(uint8(255));   
% [图片变换]保存图像按钮
function save1_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
if FileName==0
    return;
else
    h=getframe(handles.pic2);
    imwrite(h.cdata,[PathName,FileName]);
end
% [图片处理]清除图片按钮
function clear2_Callback(hObject, eventdata, handles)
    global val
    val=0;
axes(handles.pic3);
imshow(uint8(255));
% [图片处理]保存图像按钮
function save2_Callback(hObject, eventdata, handles)
    [FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
if FileName==0
    return;
else
    h=getframe(handles.pic3);
    imwrite(h.cdata,[PathName,FileName]);
end
