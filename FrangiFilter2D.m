function [outIm,whatScale,OI] = FrangiFilter2D(I, Options)
% % This function FRANGIFILTER2D uses the eigenvectors of the Hessian to
% % compute the likeliness of an image region to vessels, according
% % to the method described by Frangi:2001 (Chapter 2).
% %
% % [J,Scale,Direction] = FrangiFilter2D(I, Options)
% %
% % inputs,
% %   I : The input image (vessel image)
% %   Options : Struct with input options,
% %       .FrangiScaleRange : The range of sigmas used, default [1 10]尺度范围,最小值决定能的感应最下血管范围,最大值决定能的感应最大血管范围
% %       .FrangiScaleRatio : Step size between sigmas, default 2  每个尺度
% %       .FrangiBetaOne : Frangi correction constant, default 0.5
% %       .FrangiBetaTwo : Frangi correction constant, default 15 越小,对灰度的感应越灵敏,噪点越多
% %       .BlackWhite : Detect black ridges (default) set to true, for
% %                       white ridges set to false.
% %       .verbose : Show debug information, default true
% %
% % outputs,
% %   J : The vessel enhanced image (pixel is the maximum found in all scales)
% %   whatScale : Matrix with the scales on which the maximum intensity 
% %           of every pixel is found
% %   whatScale :取出的值中,是哪个尺度下的值
% %   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)   
% %
% % Example,
% %   I=double(imread ('vessel.png'));
% %   Ivessel=FrangiFilter2D(I);
% %   figure,
% %   subplot(1,2,1), imshow(I,[]);
% %   subplot(1,2,2), imshow(Ivessel,[0 0.25]);
% %
% % Written by Marc Schrijver, 2/11/2001
% % Re-Written by D.Kroon University of Twente (May 2009)

%clc,clear,close all
Img=I;
%Img=imread ('08_test.tif');
I=double(Img(:,:,1));%黑色管状目标
I=padarray(I,[10,10],'replicate','both');
defaultoptions = struct('FrangiScaleRange', [1 10], 'FrangiScaleRatio', 2, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo',10, 'verbose',true,'BlackWhite',true);

% Process inputs
if(~exist('options','var')) 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options)))
        warning('FrangiFilter2D:unknownoption','unknown options found');
    end
end

sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);%从[1 10]步长为2
sigmas = sort(sigmas, 'ascend');%从小到大排序

beta  = 2*options.FrangiBetaOne^2;
c     = 2*options.FrangiBetaTwo^2;

% Make matrices to store all filterd images
[nx,ny]=size(I);
ALLfiltered=zeros([size(I) length(sigmas)]);%产生sigmas个I大小的矩阵
ALLangles=zeros([size(I) length(sigmas)]);
% Frangi filter for all sigmas
for i = 1:length(sigmas)
    % Show progress
%     if(options.verbose)
%         disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]);
%     end
    
    % Make 2D hessian
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));
    
    % Correct for scale
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
   
    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);%求特征值Lambda2(较小),Lambda1,和向量
    % Compute the direction of the minor eigenvector计算较小特征向量方向
    angles = atan2(Ix,Iy);

    % Compute some similarity measures相似性侧度计算
    Lambda1(Lambda1==0) = eps;%分母不能为0
    Rb = (Lambda2./Lambda1).^2;       %当为管状结构时，Rb趋于0，点状结构时，趋于1(平滑所有有边缘效应)
    S2 =Lambda1.^2 + Lambda2.^2;  %当为噪声时，S2很大exp(-S2)趋于0 
   
    % Compute the output image
    Ifiltered = exp(-Rb/beta) .*(1-exp(-S2/c));
    % see pp. 45
    if(options.BlackWhite)
        Ifiltered(Lambda1<0)=0;%如果背景是黑的,就取大特征值小于零的,其余变为0为背景
    else
        Ifiltered(Lambda1>0)=0;%如果目标是黑的,就取大特征值大于零的,小于零的地方变为0为背景
    end
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
end

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1
    [outIm,whatScale] = max(ALLfiltered,[],3);%取各个尺度中对应位置最大值
    
%     outIm =uint8(max(max(I))*(outIm>0));
%     outIm =uint8(max(max(I))*(outIm>0));
%     outIm =uint8(max(max(I))*(outIm>max(max(outIm))/3 ));
     
    outIm = reshape(outIm,size(I));
%     if(nargout>1)
%         whatScale = reshape(whatScale,size(I));
%     end
%     if(nargout>2)
%         Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
%     end
else
    outIm = reshape(ALLfiltered,size(I));
%     if(nargout>1)
%             whatScale = ones(size(I));
%     end
%     if(nargout>2)
%         Direction = reshape(ALLangles,size(I));
%     end
end
OI=outIm(11:nx-10,11:ny-10);
%figure, imshow(uint8(Img));title('原图')
%figure, imshow(OI./max(max(OI)));title('结果')
