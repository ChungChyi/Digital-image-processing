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
axes(handles.pic1);                                        %��ʼ��ԭʼͼ��
imshow(uint8(255));
axes(handles.pic2);                                    %��ʼ������ͼ��
imshow(255*ones(255,255));
axes(handles.pic3);                                       %��ʼ������ͼ��
imshow(255*ones(255,255));
set(handles.NoiseChoose,'Enable','off');                          %ͼƬδѡ��ʱ���˲��������������ѡ��Աȷ�������ѡ
set(handles.CoChange,'Enable','off');
set(handles.HisSta,'Enable','off');
set(handles.FreFiltering,'Enable','off');
set(handles.SpaFiltering,'Enable','off');
set(handles.ColorEn,'Enable','off');
global Size ProNimg Proimg img val;
Size=3;ProNimg=[255*ones(255,255)];            %ȫ�ֱ����ĳ�ʼ����SizeΪģ���С��ProNimgΪ���������ͼ��
img=[255*ones(255,255)];                                 %imgΪѡ��ͼƬ��ProimgΪ�˲���ͼ��
Proimg=[255*ones(255,255)];
val=0;   

function varargout = DGP_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
% --------------------------------------------------------------------
%%�ļ��˵�
function file_Callback(hObject, eventdata, handles)
%�ļ�-�򿪺���
function Open_Callback(hObject, eventdata, handles)
global img;
[name,path]=uigetfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                'ѡ��ͼƬ');                    %��ȡͼ��
img=imread([path,name]);
axes(handles.pic1);
imshow(img);                                                      %��ʾͼ�񲢿�������ѡ�˵�
set(handles.NoiseChoose,'Enable','on');                          %ͼƬδѡ��ʱ���˲��������������ѡ��Աȷ�������ѡ
set(handles.CoChange,'Enable','on');
set(handles.HisSta,'Enable','on');
set(handles.FreFiltering,'Enable','on');
set(handles.SpaFiltering,'Enable','on');
set(handles.ColorEn,'Enable','on');
%�ļ�-�˳�����
function Quit_Callback(hObject, eventdata, handles)
close(gcf);

% --------------------------------------------------------------------
%%����ѡ��˵�
function NoiseChoose_Callback(hObject, eventdata, handles)   
%%����ѡ��-��˹���� �˵�
function GN_Callback(hObject, eventdata, handles)
   %%Ĭ�Ϸ���0.01 ����
   function GMFC_Callback(hObject, eventdata, handles)
global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'gaussian',0,0.01));
axes(handles.pic2);
imshow(ProNimg);   
   %%��ѡ���� ����
   function GAC_Callback(hObject, eventdata, handles)
     global img ProNimg val;
val=1;
str=['�������˹��������'];
title=['��˹������ѡ����'];
input=inputdlg(str,title,1,{'0.01'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'gaussian',0,input));
axes(handles.pic2);
imshow(ProNimg);
%%����ѡ��-������� �˵�
function RN_Callback(hObject, eventdata, handles)
    %%Ĭ�Ϸ���0.01 ����
    function RMFC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'speckle',0.01));
axes(handles.pic2);
imshow(ProNimg);    
    %%��ѡ����  ����
    function RAC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
str=['�����������������'];
title=['�������������ѡ����'];
input=inputdlg(str,title,1,{'0.04'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'speckle',input));
axes(handles.pic2);
imshow(ProNimg);
%%����ѡ��-�������� �˵�
function SN_Callback(hObject, eventdata, handles)
    %%Ĭ��ǿ��0.05 ����
    function MS_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'salt & pepper',0.05));
axes(handles.pic2);
imshow(ProNimg);
    %%��ѡ���� ����
    function SAC_Callback(hObject, eventdata, handles)
    global img ProNimg val;
val=1;
str=['�����뽷������ǿ��'];
title=['����������ѡ����'];
input=inputdlg(str,title,1,{'0.05'});
input=cell2mat(input);
input=str2num(input);
ProNimg=uint8(imnoise(img,'salt & pepper',input));
axes(handles.pic2);
imshow(ProNimg);    
%%����ѡ��-�������� ����
    function PN_Callback(hObject, eventdata, handles)
global img ProNimg val;
val=1;
ProNimg=uint8(imnoise(img,'poisson'));
axes(handles.pic2);
imshow(ProNimg);    



% --------------------------------------------------------------------
%%����任�˵�
function CoChange_Callback(hObject, eventdata, handles)
  %%����任-�Ҷȱ任 ����
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
%%ֱ��ͼͳ�Ʋ˵�
function HisSta_Callback(hObject, eventdata, handles)
    %%ֱ��ͼͳ-Rֱ��ͼ ����
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
    %%ֱ��ͼͳ-Gֱ��ͼ ����
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
    %%ֱ��ͼͳ-Bֱ��ͼ ����
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
    %%ֱ��ͼͳ-ֱ��ͼ ����
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
%%Ƶ�����˲���ǿ
function FreFiltering_Callback(hObject, eventdata, handles)
    %%Ƶ���˲�-Ƶ��ͼ ����
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
    %%Ƶ���˲�-������˹��ͨ�˲� ����
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
    %%Ƶ���˲�-������˹��ͨ�˲� ����
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
    %%Ƶ���˲�-��ɢ����Ҷ�任 ����
    function DFF_Callback(hObject, eventdata, handles)
        global img 
    axes(handles.pic3);
    i=im2double(img);
    f1=fft2(i);
    fc1=fftshift(f1);
    i=log(1+abs(fc1));
    imshow(i);
    xlabel('The Fourier transform');
    %%Ƶ���˲�-��ɢ���ұ任 ����
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
%%�ռ����˲���ǿ �˵�
function SpaFiltering_Callback(hObject, eventdata, handles)
%%�ռ����˲�-ƽ���˲����� �˵�
function PF_Callback(hObject, eventdata, handles)
    %%�ռ����˲�-ƽ���˲�����-��ֵ�˲� ����
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
    %%�ռ����˲�-ƽ���˲�����-��ֵ�˲� ����
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
    %%�ռ����˲�-ƽ���˲�����-����Ӧ�˲� ����
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
%%�ռ����˲�-���˲� �˵�
function RF_Callback(hObject, eventdata, handles)
    %%�ռ����˲�-���˲�-�ݶ��˲� ����
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
    %%�ռ����˲�-���˲�-�ݶ��� ����
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
    %%�ռ����˲�-���˲�-������˹�����˲� ����
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
    %%�ռ����˲�-���˲�-������˹������ ���� 
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
    %%�ռ����˲�-���˲�-�����˲� �˵�
    function DF_Callback(hObject, eventdata, handles)
        %%�ռ����˲�-���˲�-�����˲�-��ֱ�����˲� ����
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
        %%�ռ����˲�-���˲�-�����˲�-ˮƽ�����˲� ����
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
        %%�ռ����˲�-���˲�-�����˲�-�ԽǶ����˲� ����    
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
%%��ɫ��ǿ�˵�
function ColorEn_Callback(hObject, eventdata, handles)
    %%��ɫ��ǿ-���ɫ��ǿ �˵�
    function TCH_Callback(hObject, eventdata, handles)
        %%��ɫ��ǿ-���ɫ��ǿ-��ɫ���� ����
        function TRed_Callback(hObject, eventdata, handles)
        global img   
axes(handles.pic3);
imshow(img(:,:,1));
xlabel('red component');
        %%��ɫ��ǿ-���ɫ��ǿ-��ɫ���� ����    
        function TBlue_Callback(hObject, eventdata, handles)
         global img       
axes(handles.pic3);
imshow(img(:,:,2));
xlabel('blue component');
        %%��ɫ��ǿ-���ɫ��ǿ-��ɫ���� ����
        function TGreen_Callback(hObject, eventdata, handles)
                   global img       
axes(handles.pic3);
imshow(img(:,:,3));
xlabel('green component');
%%��ɫ��ǿ-α��ɫ��ǿ �˵�
function FCH_Callback(hObject, eventdata, handles)
    %%��ɫ��ǿ-α��ɫ��-�ܶȷָ� ����
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
    %%��ɫ��ǿ-α��ɫ��-�Ҷȼ� ����
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
    
    
% ѡ��ͼƬ��ť
function choose_Callback(hObject, eventdata, handles)
    global img;
[name,path]=uigetfile({'*.jpg';'*.bmp'},'ѡ��ͼƬ');                    %��ȡͼ��
img=imread([path,name]);
axes(handles.pic1);
imshow(img);                                                      %��ʾͼ�񲢿�������ѡ�˵�
set(handles.NoiseChoose,'Enable','on');                          %ͼƬδѡ��ʱ���˲��������������ѡ��Աȷ�������ѡ
set(handles.CoChange,'Enable','on');
set(handles.HisSta,'Enable','on');
set(handles.FreFiltering,'Enable','on');
set(handles.SpaFiltering,'Enable','on');
set(handles.ColorEn,'Enable','on');
%[ͼƬ�任]���ͼƬ��ť
function clear1_Callback(hObject, eventdata, handles)
    global val
    val=0;
axes(handles.pic2);
imshow(uint8(255));   
% [ͼƬ�任]����ͼ��ť
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
% [ͼƬ����]���ͼƬ��ť
function clear2_Callback(hObject, eventdata, handles)
    global val
    val=0;
axes(handles.pic3);
imshow(uint8(255));
% [ͼƬ����]����ͼ��ť
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
