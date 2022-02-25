function dfundyNum = Dy(img,npx,npy,dy,n,ooa)

Dy = FDMatrixY(npx,npy,dy,n,ooa);
dfundyNum = Dy*img(:);
dfundyNum = reshape(dfundyNum,npx,npy);

end