% clc
% clear all
% close all

load('Brain2D.mat');
FOV=256;
Nc = 12;
Nx =  FOV;
Ny =  FOV;
rate =2;


 M=1;
% mask=zeros(Nx,Ny);
% mask(1:M:end,:)=1;
mask=k;

%% Images of each coils

cim=ifftshift(ifft2(ifftshift(DATA)));


c_raw=DATA;

squared_img = power(abs(cim), 2);
sum_img = sum(squared_img, 3);
rsos = sqrt(sum_img);
figure,
imshow(abs(rsos),[])
%% different noise

min_a = min(min(c_raw(:)));
max_a = max(max(c_raw(:)));
for n=1:Nc
    normal(:,:,n) = (c_raw(:,:,n)-min_a)./abs(max_a-min_a); 
end 

%% Coil images
coil_image=ifftshift(ifft2(ifftshift(normal)));


figure(1),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(coil_image(:,:,n)),[])
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
 

%% Coil images for coil sensitivity from sample (filtered)
for n=1:Nc
    coil_image1(:,:,n)=ifftshift(ifft2(ifftshift(padded_filter(:,:,n))));   
end
figure(1),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(coil_image1(:,:,n)),[])
end

%% Combined SOS image (low resolution)for coil sensitivty
for n=1:Nc
squared_img(:,:,n) = power(abs(coil_image1(:,:,n)), 2);
end

sum_img = sum(squared_img, 3);
rsos1 = sqrt(sum_img);
figure,
imshow((abs(rsos1)),[])

%% Reference Image for error calculation and difference Image
for n=1:Nc
sq_img(:,:,n) = power(abs(coil_image(:,:,n)), 2);
end
s_img = sum(sq_img, 3);
xo = sqrt(s_img);
figure,
imshow((abs(xo)),[])

%% Coil sensitivity

for n=1:Nc
c_sens2(:,:,n)=coil_image1(:,:,n)./rsos1;
end
figure,
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(c_sens2(:,:,n)),[])
end




normal = mask.*normal;
%% ifft noisy k space
image_domain_with_noise=zeros(Nx,Ny,Nc);
image_domain_with_noise = (ifftshift(ifft2(ifftshift(normal)))); 
figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(image_domain_with_noise(:,:,n)),[])
end







conjugate_sens=conj(c_sens2);

sum_for_a=zeros(Nx,Ny,Nc);
for n=1:Nc
sum_for_a(:,:,n)=conjugate_sens(:,:,n).*image_domain_with_noise(:,:,n);
end
figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(sum_for_a(:,:,n)),[])
end





a=zeros(Nx,Ny);
    for n=1:Nc
        a=a+sum_for_a(:,:,n);
    end
figure,
imshow(abs(a),[])




b=zeros(Nx,Ny,101);
b(:,:,1)=0;
p=a;
r=zeros(Nx,Ny,101);
r(:,:,1)=a;


for i=1:100
    for n=1:Nc
       Ewithp(:,:,n)=c_sens2(:,:,n).*p;
    end
    k_space=fftshift(fft2(fftshift(Ewithp)));
      k_sample=mask.*k_space;
    Eh_half=ifftshift(ifft2(ifftshift(k_sample)));
    Eh=conjugate_sens.*Eh_half;
      q=Eh;
      qsum=zeros(Nx,Ny);
    for n=1:Nc
        qsum=qsum+q(:,:,n);
    end
permute_r=permute(r,[1 2 3]);
reshape_r=reshape(permute_r,256*256,[]);
    alpha=reshape_r(:,i)'*reshape_r(:,i)./(p(:)'*qsum(:));
    b(:,:,i+1)=b(:,:,i)+alpha*p;
    r(:,:,i+1)=r(:,:,i)-alpha*qsum;
    beta=reshape_r(:,i+1)'*reshape_r(:,i+1)./(reshape_r(:,i)'*reshape_r(:,i));
    p=r(:,:,i+1)+beta*p;
 end

figure,
imshow(abs(b(:,:,100)),[])
j=[2:101];
delta=sum(sum(abs(b(:,:,j)-rsos)./(abs(rsos)+eps)));
sca=zeros(1,100);
sca(:,:)=delta(:,:,:);
figure,
plot(j,sca)
% hold on
% j=[2:101];
% delta=sum(sum(abs(b(:,:,j)-rsos)./(abs(rsos)+eps)));
% scal=zeros(1,100);
% scal(:,:)=delta(:,:,:);
% figure,
% plot(j,sca,'--o')
% hold on
% plot(j,scal,'--*')
% hold on
% j=[2:101];
% delta=sum(sum(abs(b(:,:,j)-rsos)./(abs(rsos)+eps)));
% scala=zeros(1,100);
% scala(:,:)=delta(:,:,:);
% plot(j,scala,'--^')
% legend('M=8','M=4','M=2')


diff=abs(xo)-abs(b(:,:,21));
figure,
imshow(abs(diff),[])
imagesc(abs(diff))

error2 = (abs(xo)-abs(b(:,:,10))).^2;
RMSE2 = sqrt(sum(error2(:))/(Nx * Ny));
NRMSE2 = RMSE2/(Nx*Ny)
% 
% figure,
% imshow(abs(qsum),[])