
%% phantom data
load 'modifiedshep.mat'
FOV=256;
ph=phantom('modified shepp-logan',FOV);
figure(1) ;
imshow(ph,[])

%%

IMSIZE=256*256;

img_WT=FWT2_PO(xo_u,6,MakeONFilter('Daubechies',12));
figure,
imshow(abs(img_WT),[])

index_WT=sort(abs(img_WT(:)),1,'descend');
c=0;

for pctg=floor([0.01,0.05,0.1,0.2,0.3,0.5]*IMSIZE);
    c=c+1;
    thresh(:,c)=index_WT(pctg);
    tmp=img_WT;
    tmp(find(abs(tmp)<thresh(:,c)))=0;
    rec_wav(:,:,c)=IWT2_PO(tmp,6,MakeONFilter('Daubechies',12));
    
end

tmp = [];
for n=1:6
tmp = cat(2,tmp,abs(cat(1, rec_wav(:,:,n))));
end
figure, imshow(tmp,[],'InitialMagnification',100), title('1% , 5% , 10% , 20% , 30% , 50% '), 
ylabel(' Wavelet'), drawnow,

disp('Done')


