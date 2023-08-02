clc
clear all
close all

addpath(strcat(pwd,'/utils'));

if exist('FWT2_PO') <2
	error('must have Wavelab installed and in the path');
end
%% Load file
load('modifiedshep.mat');

%% Parameters
FOV=256;
ph=phantom('modified shepp-logan',FOV);
Nc = 8;
Nx =  FOV;
Ny =  FOV;

%% Images of each coils
% 
for n=1:Nc
    c_img1(:,:,n) = ph.*c_sens(:,:,n); 
end

c_raw2=fftshift(fft2(fftshift(c_img1)));

squared_img = power(abs(c_img1), 2);
sum_img = sum(squared_img, 3);
sos = sqrt(sum_img);
figure,
imshow(abs(sos),[])

normal=c_raw2;
%% Normalization
% min_a = min(min(DATA(:)));
% max_a = max(max(DATA(:)));
% for n=1:Nc
%     normal(:,:,n) = (DATA(:,:,n)-min_a)./abs(max_a-min_a); 
% end 
% 
% min_a = min(min(normal(:)));
% max_a = max(max(normal(:)));
% for n=1:Nc
%     normal(:,:,n) = (normal(:,:,n)-min_a)./abs(max_a-min_a); 
% end 

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
coil_sens(:,:,n)=coil_img1(:,:,n)./rsos1;
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
TVWeight = 3.5; 	% Weight for TV penalty
xfmWeight =0;	% Weight for Transform L1 penalty
Itnlim = 5;	
pctg = [0.5];  	% undersampling factor
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
              B(L,1:M)=coil_sens(y:delta:end,x,L);
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
b=zeros(Nx,Ny,50);

for n=1:50
	res = mytryfnlCgPI(res,param,coil_sens,normal,mask,Nc);
	im_res = XFM'*res;
figure(10), imshow(abs(im_res),[])
 b(:,:,n)=im_res;
end


% figure, imshow(abs(cat(2,im_dc,im_res)),[]);
% figure, imshow(abs(cat(2,im_dc(155:304,110:209), im_res(155:304,110:209))),[0,1],'InitialMagnification',200);
% title(' zf-w/dc              l_1 Wavelet');

%%

% diff_image=sos-im_res;
% figure,
% imagesc(abs(diff_image))

error = (abs(sos)-abs(im_res)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

% j=[2:100];
% delta=sum(sum(abs(b(:,:,j)-sos)./(abs(sos)+eps)));
% sca=zeros(1,99);
% sca(:,:)=delta(:,:,:);
% figure,
% plot(j,sca)


% j=[2:100];
% delta=sum(sum(abs(b(:,:,j)-image)./(abs(image)+eps)));
% scal=zeros(1,99);
% scal(:,:)=delta(:,:,:);
% figure,
% plot(j,scal)
% hold on
% plot(j,sca)
% legend('M=4','M=8')