N = [Nx,Ny];	% image Size
DN = [Nx,Ny];	% data Size	% Weight for Transform L1 penalty
Itnlim = 5;	
pctg = [0.5];  	% undersampling factor
P = 5;%%5
M=1;

pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
mask = genSampling(pdf,10,60);		% generates a sampling pattern

normal=(4.*normal).*mask;

k=1;
l=1;

for j=[0 1e-2 1.1e-3 1.5e-4 2e-5 1e-6 1e-7 2.5e-7 1e-8]
for i=[0 1e-2 1.1e-3 1.5e-4 2e-5 1e-6 1e-7 2.5e-7 1e-8]
TVWeight = i%1e-7; 	% Weight for TV penalty
xfmWeight = j%2.5e-7

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

% figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
for n=1:100
	res = mytryfnlCgPI(res,param,c_sens,normal,mask,Nc);
	im_res = XFM'*res;
% figure(10), imshow(abs(im_res),[])
end

res=1e-3.*abs(im_res);

error = ((sos)-abs(res)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

func(k,l)=NRMSE;
k=k+1;
end
k=1;
l=l+1;
end
 
imagesc(func)
colorbar
set(gca,'YDir','normal');
% xlim([0 0.02])
% ylim([0 0.02])

 