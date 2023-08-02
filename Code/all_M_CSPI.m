clc
clear all
close all

addpath(strcat(pwd,'/utils'));

if exist('FWT2_PO') <2
	error('must have Wavelab installed and in the path');
end
%% Load file
load('Brain2D');

%% Parameters
FOV=256;
Nc = 12;
Nx =  FOV;
Ny =  FOV;

%% Normalization
min_a = min(min(DATA(:)));
max_a = max(max(DATA(:)));
for n=1:Nc
    normal(:,:,n) = (DATA(:,:,n)-min_a)./abs(max_a-min_a); 
end 

%% Coil images
coil_img=ifftshift(ifft2(ifftshift(normal)));


% figure(1),
% for n=1:Nc
%     subplot(2,ceil(Nc/2),n)
%     imshow(abs(coil_img(:,:,n)),[])
% end

%% 32x32 sample

Sample_from_k=zeros(40,40,Nc);
Sample_from_k(:,:,:)=normal(109:148,109:148,:);
 
% Hamming window filteration
w = hann(40)*hann(40)';
FFT = fftshift(fft2(w)); % complex matrix
FFT_abs = abs(FFT); % absolut values of fft
imagesc(1+log10(FFT_abs)) % show fft results
w_new = ifft2(ifftshift(FFT)); % you need fft not fft_abs
% figure,
% imshow(abs(w_new),[])


% Dot multiplication of window with sample

for n=1:Nc
filtered_sample(:,:,n)=w_new.*Sample_from_k(:,:,n);
end

% fitting 32x32 in back and zero padding
padded_filter=zeros(Nx,Ny,Nc);
padded_filter(109:148,109:148,:)=filtered_sample(:,:,:);
 

% Coil images for coil sensitivity from sample (filtered)
for n=1:Nc
    coil_img1(:,:,n)=ifftshift(ifft2(ifftshift(padded_filter(:,:,n))));   
end
% figure(1),
% for n=1:Nc
%     subplot(2,ceil(Nc/2),n)
%     imshow(abs(coil_img1(:,:,n)),[])
% end

% Combined SOS image (low resolution)for coil sensitivty
for n=1:Nc
squared_img(:,:,n) = power(abs(coil_img1(:,:,n)), 2);
end

sum_img = sum(squared_img, 3);
rsos1 = sqrt(sum_img);
% figure,
% imshow((abs(rsos1)),[])


%% coil sensitivity
for n=1:Nc
c_sens(:,:,n)=coil_img1(:,:,n)./rsos1;
end
% figure,
% for n=1:Nc
%     subplot(2,ceil(Nc/2),n)
%     imshow(abs(c_sens(:,:,n)),[])
% end



%% Reference Image for error calculation and difference Image
for n=1:Nc
sq_img(:,:,n) = power(abs(coil_img(:,:,n)), 2);
end
s_img = sum(sq_img, 3);
sos = sqrt(s_img);
figure,
imshow((abs(sos)),[])


 N = [Nx,Ny];	% image Size
DN = [Nx,Ny];	% data Size
TVWeight = 1e-7; 	% Weight for TV penalty
xfmWeight =2.5e-7;	% Weight for Transform L1 penalty
Itnlim = 5;	
pctg = [0.5];  	% undersampling factor
P = 5;%%5
M=1;

pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
mask = genSampling(pdf,10,60);		% generates a sampling pattern

normal=(2.*normal).*mask;

noise_along_y=zeros(Nx,Ny,Nc);
for n=1:Nc
noise_along_y(:,:,n) = (ifftshift(ifft2(ifftshift(normal(:,:,n))))); 
end
figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(noise_along_y(:,:,n)),[])
end


 delta=Ny/M;

recon_imge=zeros(Nx,Ny);

for x=1:Nx
     for y=1:delta
          for L=1:Nc
              B(L,1:M)=c_sens(y:delta:end,x,L);
              pixel_vector(L,1)=noise_along_y(y,x,L);
          end
          inv_B=pinv(B);
          recon_imge(y:delta:end,x)=inv_B*pixel_vector;
     end
end

figure,
imshow(abs(recon_imge),[])

image=recon_imge;
%%
FT = p2DFT(mask, N, 1, 2);
data = FT*image;

im_dc = FT'*(data.*mask./pdf);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
b=zeros(Nx,Ny,150);

for n=1:100
	res = mytryfnlCgPI(res,param,c_sens,normal,mask);
	im_res = XFM'*res;
figure(10), imshow(abs(im_res),[])
 b(:,:,n)=im_res;
end


figure, imshow(abs(cat(2,im_dc,im_res)),[]);

error = (abs(sos)-abs(im_res)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

j=[2:100];
delta=sum(sum(abs(b(:,:,j)-sos)./(abs(sos)+eps)));
sca=zeros(1,99);
sca(:,:)=delta(:,:,:);

%%
 N = [Nx,Ny];	% image Size
DN = [Nx,Ny];	% data Size
TVWeight = 2e-6; 	% Weight for TV penalty
xfmWeight =5e-7;	% Weight for Transform L1 penalty
Itnlim = 5;	
pctg = [0.25];  	% undersampling factor
P = 5;%%5
M=1;

pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
mask = genSampling(pdf,10,60);		% generates a sampling pattern

normal=(4.*normal).*mask;

noise_along_y=zeros(Nx,Ny,Nc);
for n=1:Nc
noise_along_y(:,:,n) = (ifftshift(ifft2(ifftshift(normal(:,:,n))))); 
end
figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(noise_along_y(:,:,n)),[])
end


 delta=Ny/M;

recon_imge=zeros(Nx,Ny);

for x=1:Nx
     for y=1:delta
          for L=1:Nc
              B(L,1:M)=c_sens(y:delta:end,x,L);
              pixel_vector(L,1)=noise_along_y(y,x,L);
          end
          inv_B=pinv(B);
          recon_imge(y:delta:end,x)=inv_B*pixel_vector;
     end
end

figure,
imshow(abs(recon_imge),[])

image=recon_imge;
%%

FT = p2DFT(mask, N, 1, 2);
data = FT*image;

im_dc = FT'*(data.*mask./pdf);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
b=zeros(Nx,Ny,150);

for n=1:100
	res = mytryfnlCgPI(res,param,c_sens,normal,mask);
	im_res = XFM'*res;
figure(10), imshow(abs(im_res),[])
 b(:,:,n)=im_res;
end


figure, imshow(abs(cat(2,im_dc,im_res)),[]);

error = (abs(sos)-abs(im_res)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

j=[2:100];
delta=sum(sum(abs(b(:,:,j)-sos)./(abs(sos)+eps)));
scal=zeros(1,99);
scal(:,:)=delta(:,:,:);



%%
 N = [Nx,Ny];	% image Size
DN = [Nx,Ny];	% data Size
TVWeight = 4e-6; 	% Weight for TV penalty
xfmWeight =10e-7;	% Weight for Transform L1 penalty
Itnlim = 5;	
pctg = [0.125];  	% undersampling factor
P = 5;%%5
M=1;

pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
mask = genSampling(pdf,10,60);		% generates a sampling pattern

normal=(8.*normal).*mask;

noise_along_y=zeros(Nx,Ny,Nc);
for n=1:Nc
noise_along_y(:,:,n) = (ifftshift(ifft2(ifftshift(normal(:,:,n))))); 
end
figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(noise_along_y(:,:,n)),[])
end


 delta=Ny/M;

recon_imge=zeros(Nx,Ny);

for x=1:Nx
     for y=1:delta
          for L=1:Nc
              B(L,1:M)=c_sens(y:delta:end,x,L);
              pixel_vector(L,1)=noise_along_y(y,x,L);
          end
          inv_B=pinv(B);
          recon_imge(y:delta:end,x)=inv_B*pixel_vector;
     end
end

figure,
imshow(abs(recon_imge),[])

image=recon_imge;
%%

FT = p2DFT(mask, N, 1, 2);
data = FT*image;

im_dc = FT'*(data.*mask./pdf);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
b=zeros(Nx,Ny,150);

for n=1:100
	res = mytryfnlCgPI(res,param,c_sens,normal,mask);
	im_res = XFM'*res;
figure(10), imshow(abs(im_res),[])
 b(:,:,n)=im_res;
end


figure, imshow(abs(cat(2,im_dc,im_res)),[]);

error = (abs(sos)-abs(im_res)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

j=[2:100];
delta=sum(sum(abs(b(:,:,j)-sos)./(abs(sos)+eps)));
scala=zeros(1,99);
scala(:,:)=delta(:,:,:);


figure,
plot(j,sca)
hold on
plot(j,scal)
hold on
plot(j,scala)
legend('2','4','8')
