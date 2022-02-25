function dfundxNum = Dx(img,npx,npy,dx,n,ooa)

Dx = FDMatrixX(npx,npy,dx,n,ooa);
dfundxNum = Dx*img(:);
dfundxNum = reshape(dfundxNum,npx,npy);

end

