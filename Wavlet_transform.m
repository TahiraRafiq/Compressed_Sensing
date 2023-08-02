B=xo;
wavelet='db1';
level=2;
n=2;
rv=256;
w='db1';
[C,S]=wavedec2(B,n,w);

A = cell(1,level); H = A; V = A; D = A;

for k = 1:level
    A{k} = appcoef2(C,S,wavelet,k); % approx
    [H{k} V{k} D{k}] = detcoef2('a',C,S,k); % details  
    A{k} = wcodemat(A{k},rv);
    H{k} = wcodemat(H{k},rv);
    V{k} = wcodemat(V{k},rv);
    D{k} = wcodemat(D{k},rv);
end
    dec = cell(1,level);
    dec{level} = [A{level} H{level} ; V{level} D{level}];
    
    for k = level-1:-1:1
        dec{k} = [imresize(dec{k+1},size(H{k})) H{k} ; V{k} D{k}];
    end
    
    image(dec{1});
    colormap((gray))
    


