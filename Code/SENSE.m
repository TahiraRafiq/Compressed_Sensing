% clc
% clear all
% close all

%% Load file
load('Brain2D');

%% Parameters
FOV=256;
Nc = 12;
Nx =  FOV;
Ny =  FOV;
rate = 2;
M=1;
% mask=zeros(Nx,Ny);
% mask(1:M:end,:)=1;
% mask=k;

min_a = min(min(DATA(:)));
max_a = max(max(DATA(:)));
for n=1:Nc
    normal(:,:,n) = (DATA(:,:,n)-min_a)./(max_a-min_a); 
end
% norm=DATA;

%% Coil images
coil_img=ifftshift(ifft2(ifftshift(normal)));
figure(1),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(coil_img(:,:,n)),[])
end

%% 32x32 sample

Sample_from_k=zeros(40,40,Nc);
Sample_from_k(:,:,:)=normal(109:148,109:148,:);
 
%% Hamming window filteration
w = hann(40)*hann(40)';
FFT = fftshift(fft2(w)); % complex matrix
FFT_abs = abs(FFT); % absolut values of fft
imagesc(1+log10(FFT_abs)) % show fft results
w_new = ifft2(ifftshift(FFT)); % you need fft not fft_abs
figure,
imshow(abs(w_new),[])


% Dot multiplication of window with sample

for n=1:Nc
filtered_sample(:,:,n)=w_new.*Sample_from_k(:,:,n);
end

% fitting 32x32 in back and zero padding
padded_filter=zeros(Nx,Ny,Nc);
padded_filter(109:148,109:148,:)=filtered_sample(:,:,:);
 

%% Coil images
for n=1:Nc
    coil_img1(:,:,n)=ifftshift(ifft2(ifftshift(padded_filter(:,:,n))));   
end
figure(1),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(coil_img1(:,:,n)),[])
end

%% Combined SOS image
 
for n=1:Nc
squared_img(:,:,n) = power(abs(coil_img1(:,:,n)), 2);
end

sum_img = sum(squared_img, 3);
rsos1 = sqrt(sum_img);
figure,
imshow((abs(rsos1)),[])


for n=1:Nc
sq_img(:,:,n) =power(abs(coil_img(:,:,n)),2);
end
s_img = sum(sq_img, 3);
xo = sqrt(s_img);
figure,
imshow((abs(xo)),[])


%% Coil sensitivity

for n=1:Nc
c_sens(:,:,n)=coil_img1(:,:,n)./rsos1;
end
figure,
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(c_sens(:,:,n)),[])
end


normal = mask.*normal;

%% ifft noisy k space

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
              pixel_vector(L,1)=noise_along_y(y,x,L)*M;
          end
          inv_B=pinv(B);
          recon_imge(y:delta:end,x)=inv_B*pixel_vector;
     end
end

figure,
imshow(abs(recon_imge),[])




error = (abs(xo)-abs(recon_imge)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)


diff=abs(xo)-abs(recon_imge);
figure,
imshow(abs(diff),[])


% figure,
% imshow(abs(abs(xo)-abs(recon_imge)*2),[0 2e-7])


% error2 = (xo-abs(recon_imge) ).^2;
% RMSE = sqrt(sum(sum(error))/(Nx * Ny));
% NRMSE = RMSE/(Nx*Ny);