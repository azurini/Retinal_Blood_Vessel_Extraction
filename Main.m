%Author : ����������
%
clc;
clear;
close all;
%% step1��A read ����Ĥͼ
%
Img = imread('08_test.tif');
mash = mask('08_test_mask.gif');  %��ȡ�߿򣬲�����������Ե�Թ�Ӱ��

dim = ndims(Img);
if(dim == 3)  
    inImg = rgb2gray(Img);%Input is a color image ��ͼ��ֵ��
end

   I = im2double(inImg);
   %I = imgaussian(I,2,5);
   figure, imshow(Img);title('ԭͼ');

figure, imshow(inImg);title('step1����ֵ��ͼƬ');

%% step2��B ȥ��Ѫ��
%
blood = imread('08_manual1.gif');
blood = im2double(blood);
%figure, imshow(blood);title('step1��B���ݿ���ר���ֱ��Ѫ��λ��');

bmaster = inImg;
[m, n]= find(blood==1);  %�ҳ�Ѫ��λ��
for i=1:1:length(m)
    bmaster(m(i),n(i)) = 0;
end                      %�ڶ�ֵ��ԭͼ���������Ѫ��
%figure, imshow(bmaster);title('step2���Ӷ�ֵ��ͼA��ȥ���ݿ�Ѫ��B');

[Dxx,Dxy,Dyy] = Hessian2D(I,2);
   m = im2double(mash);
   Dxx = Dxx.*m;
   Dyy = Dyy.*m;
   
bself = inImg;
[m, n]= find(Dxx>=0);  %�ҳ�Ѫ��λ��
for i=1:1:length(m)
    bself(m(i),n(i)) = 0;
end                      %�ڶ�ֵ��ԭͼ���������Ѫ��
[m, n]= find(Dyy>=0);  %�ҳ�Ѫ��λ��
for i=1:1:length(m)
    bself(m(i),n(i)) = 0;
end                      %�ڶ�ֵ��ԭͼ���������Ѫ��
figure, imshow(bself);title('step2���Ӷ�ֵ��ͼA��ȥ�Լ���⵽��Ѫ��B');

%% step3���ж��Ƿ񲡱�
%
[m, n]= find(bself>=145);    %�ҳ�����λ�ã�Ѱ�ҽ��ȡ���ڲ�������
bb=bself.*0;
for i=1:1:length(m)
    bb(m(i),n(i)) = bself(m(i),n(i));
end
bb = bb.*mash;

figure, imshow(bb);title('step3����Ǹ���λ�ã���ʱδ�������̣�');

%% step4��ȥ������

prompt1='ͨ������λ�ã��жϾ������ĸ�����1-3��\n';
a1=input(prompt1);
prompt2='�������̿�ȣ�40-100��\n';
a2=input(prompt2);

% �����Ѫ��λ��
mblood = im2double(bself);
[outIm,whatScale,OI] = FrangiFilter2D(mblood,0.5);%��������ΪѪ�ܾ�������
H = OI./max(max(OI));
figure, imshow(H);        

% �ҳ�����λ�ã�Ѱ�ҽ��ȡ���ڲ�������
[m, n]= find(H>=0.55);    
Hoptic=H.*0;              
Hdisk=[0,0];
for i=1:1:length(m)
    Hoptic(n(i),m(i)) = H(n(i),m(i));
    Hdisk(i,:) = [n(i) m(i)];
end

hold on;plot(Hdisk(:,1),Hdisk(:,2),'*');  %�궨����λ��

% ����Ѫ�ܸ�������
% k-means (���þ�������) �������ĸ����������һ���������˹��жϣ��������������������ϵľ�������
[a,b,c] = kmeans(Hdisk,a1);% bΪ������������

% ��һ���˹����� ͨ����Ѫ�ܸ�����������������
for i=1:1:a1
   hold on;plot(b(i,1),b(i,2),'o');
end
title('step4��A���Ѫ�ܸ�����������*�;�������o')

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
figure,imshow(I);title('step4��Aͨ���궨���ȥ������');

% �������Զ����� ��Ե���
BW = edge(bb,'canny');
figure, imshow(BW);title('step4��B��Ե��ⷽ��Ѱ�ҵ������̺Ͳ���');

BW2 = BW;
BW2(II<a2)=0;
figure,imshow(BW2);title('step4��Bͨ����Ե���ȥ������');


%% step5����ò���λ��
%

