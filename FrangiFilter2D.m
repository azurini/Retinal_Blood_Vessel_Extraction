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
% %       .FrangiScaleRange : The range of sigmas used, default [1 10]�߶ȷ�Χ,��Сֵ�����ܵĸ�Ӧ����Ѫ�ܷ�Χ,���ֵ�����ܵĸ�Ӧ���Ѫ�ܷ�Χ
% %       .FrangiScaleRatio : Step size between sigmas, default 2  ÿ���߶�
% %       .FrangiBetaOne : Frangi correction constant, default 0.5
% %       .FrangiBetaTwo : Frangi correction constant, default 15 ԽС,�ԻҶȵĸ�ӦԽ����,���Խ��
% %       .BlackWhite : Detect black ridges (default) set to true, for
% %                       white ridges set to false.
% %       .verbose : Show debug information, default true
% %
% % outputs,
% %   J : The vessel enhanced image (pixel is the maximum found in all scales)
% %   whatScale : Matrix with the scales on which the maximum intensity 
% %           of every pixel is found
% %   whatScale :ȡ����ֵ��,���ĸ��߶��µ�ֵ
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
I=double(Img(:,:,1));%��ɫ��״Ŀ��
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

sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);%��[1 10]����Ϊ2
sigmas = sort(sigmas, 'ascend');%��С��������

beta  = 2*options.FrangiBetaOne^2;
c     = 2*options.FrangiBetaTwo^2;

% Make matrices to store all filterd images
[nx,ny]=size(I);
ALLfiltered=zeros([size(I) length(sigmas)]);%����sigmas��I��С�ľ���
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
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);%������ֵLambda2(��С),Lambda1,������
    % Compute the direction of the minor eigenvector�����С������������
    angles = atan2(Ix,Iy);

    % Compute some similarity measures�����Բ�ȼ���
    Lambda1(Lambda1==0) = eps;%��ĸ����Ϊ0
    Rb = (Lambda2./Lambda1).^2;       %��Ϊ��״�ṹʱ��Rb����0����״�ṹʱ������1(ƽ�������б�ԵЧӦ)
    S2 =Lambda1.^2 + Lambda2.^2;  %��Ϊ����ʱ��S2�ܴ�exp(-S2)����0 
   
    % Compute the output image
    Ifiltered = exp(-Rb/beta) .*(1-exp(-S2/c));
    % see pp. 45
    if(options.BlackWhite)
        Ifiltered(Lambda1<0)=0;%��������Ǻڵ�,��ȡ������ֵС�����,�����Ϊ0Ϊ����
    else
        Ifiltered(Lambda1>0)=0;%���Ŀ���Ǻڵ�,��ȡ������ֵ�������,С����ĵط���Ϊ0Ϊ����
    end
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
end

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1
    [outIm,whatScale] = max(ALLfiltered,[],3);%ȡ�����߶��ж�Ӧλ�����ֵ
    
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
%figure, imshow(uint8(Img));title('ԭͼ')
%figure, imshow(OI./max(max(OI)));title('���')
