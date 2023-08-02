% this is a script to demonstrate the original experiment by Candes, Romberg and Tao
%
% (c) Michael Lustig 2007

rand('twister',2000);
addpath(strcat(pwd,'/utils'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [256,256]; 		% image Size
DN = [256,256]; 	% data Size
pctg = [0.5];  	% undersampling factor
P = 5;			% Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty
Itnlim = 1;		% Number of iterations


% generate variable density random sampling
pdf = genPDF(DN,P,pctg , 2 ,0.1,0); % generates the sampling PDF
%  figure,
%  imshow(abs(pdf),[]) 
k = genSampling(pdf,10,60);		% generates a sampling pattern
%  figure,
%  imshow(abs(k),[]) 


%generate image
% im = (phantom(N(1)))  + randn(N)*0.01 + i*randn(N)*0.01;

min_xo = min(min(xo(:)));
max_xo = max(max(xo(:)));
im = (xo(:,:)-min_xo)./(max_xo-min_xo); 
% figure,
%   imshow(abs(im),[]) 
%% for variable density kspace image

% image_k=fftshift(fft2(fftshift(im)));
% image_k=image_k.*k;
% image=ifftshift(ifft2(ifftshift(image_k)))
% figure,
% imshow(abs(image),[])

%% 
%generate Fourier sampling operator
FT = p2DFT(k, N, 1, 2);
data = FT*im;

%generate transform operator

XFM = Wavelet('Daubechies',6,4);	% Wavelet
%XFM = TIDCT(8,4);			% DCT
% XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

im_dc = FT'*(data./pdf);	% init with zf-w/dc (zero-fill with density compensation)
figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:5
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc


% create a low-res mask
mask_lr = genLRSampling_pctg(DN,pctg,1,0);
im_lr = ifft2c(zpad(fft2c(im).*mask_lr,N(1),N(2)));

im_full = ifft2c(zpad(fft2c(im),N(1),N(2)));
figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
title('original             low-res              zf-w/dc              TV');

figure, plot(1:N(1), abs(im_full(end/2,:)),'g',1:N(1), abs(im_lr(end/2,:)),'y', 1:N(2), abs(im_dc(end/2,:)),'m', 1:N(2), abs(im_res(end/2,:)),'b','LineWidth',2);
legend('original', 'LR', 'zf-w/dc', 'TV+WT');


% 
diff=xo-im_dc;
figure,
imagesc(abs(diff))
