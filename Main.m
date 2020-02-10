%Author : 卡拉马佐夫
%
clc;
clear;
close all;
%% step1：A read 视网膜图
%
Img = imread('08_test.tif');
mash = mask('08_test_mask.gif');  %读取边框，并用于消除边缘辉光影响

dim = ndims(Img);
if(dim == 3)  
    inImg = rgb2gray(Img);%Input is a color image 彩图二值化
end

   I = im2double(inImg);
   %I = imgaussian(I,2,5);
   figure, imshow(Img);title('原图');

figure, imshow(inImg);title('step1：二值化图片');

%% step2：B 去除血管
%
blood = imread('08_manual1.gif');
blood = im2double(blood);
%figure, imshow(blood);title('step1：B数据库中专家手标的血管位置');

bmaster = inImg;
[m, n]= find(blood==1);  %找出血管位置
for i=1:1:length(m)
    bmaster(m(i),n(i)) = 0;
end                      %在二值化原图中清除所有血管
%figure, imshow(bmaster);title('step2：从二值化图A减去数据库血管B');

[Dxx,Dxy,Dyy] = Hessian2D(I,2);
   m = im2double(mash);
   Dxx = Dxx.*m;
   Dyy = Dyy.*m;
   
bself = inImg;
[m, n]= find(Dxx>=0);  %找出血管位置
for i=1:1:length(m)
    bself(m(i),n(i)) = 0;
end                      %在二值化原图中清除所有血管
[m, n]= find(Dyy>=0);  %找出血管位置
for i=1:1:length(m)
    bself(m(i),n(i)) = 0;
end                      %在二值化原图中清除所有血管
figure, imshow(bself);title('step2：从二值化图A减去自己检测到的血管B');

%% step3：判断是否病变
%
[m, n]= find(bself>=145);    %找出高亮位置，寻找结果取决于参数设置
bb=bself.*0;
for i=1:1:length(m)
    bb(m(i),n(i)) = bself(m(i),n(i));
end
bb = bb.*mash;

figure, imshow(bb);title('step3：标记高亮位置（此时未消除视盘）');

%% step4：去除视盘

prompt1='通过高亮位置，判断聚类中心个数（1-3）\n';
a1=input(prompt1);
prompt2='设置视盘宽度（40-100）\n';
a2=input(prompt2);

% 获得主血管位置
mblood = im2double(bself);
[outIm,whatScale,OI] = FrangiFilter2D(mblood,0.5);%高亮区域为血管纠结区域
H = OI./max(max(OI));
figure, imshow(H);        

% 找出高亮位置，寻找结果取决于参数设置
[m, n]= find(H>=0.55);    
Hoptic=H.*0;              
Hdisk=[0,0];
for i=1:1:length(m)
    Hoptic(n(i),m(i)) = H(n(i),m(i));
    Hdisk(i,:) = [n(i) m(i)];
end

hold on;plot(Hdisk(:,1),Hdisk(:,2),'*');  %标定高亮位置

% 聚类血管高亮区域
% k-means (设置聚类中心) 聚类中心个数必须从上一步工作中人工判断，如果病变可能有两个以上的聚类中心
[a,b,c] = kmeans(Hdisk,a1);% b为聚类中心坐标

% （一）人工方法 通过主血管附近高亮区域标记视盘
for i=1:1:a1
   hold on;plot(b(i,1),b(i,2),'o');
end
title('step4：A获得血管附近高亮区域*和聚类中心o')

[m,n]=size(H);
II=zeros(m,n); 
for i=1:1:m
    for j=1:1:n
        II(i,j)=sqrt((i-b(1,2)).^2+(j-b(1,1)).^2);
    end
end

aa=II<a2;
I=H;
I(aa)=0;
figure,imshow(I);title('step4：A通过标定结果去除视盘');

% （二）自动方法 边缘检测
BW = edge(bb,'canny');
figure, imshow(BW);title('step4：B边缘检测方法寻找到的视盘和病变');

BW2 = BW;
BW2(II<a2)=0;
figure,imshow(BW2);title('step4：B通过边缘结果去除视盘');


%% step5：获得病变位置
%

