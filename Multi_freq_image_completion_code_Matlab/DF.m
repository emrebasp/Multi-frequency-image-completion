function dF = DF(liftedImg,theta,npx,npy, ntheta, npf,df,n,ooa)
% DEta computes the first order left-invariant derivative in \xi = \cos(theta)\partial_x + \sin(theta)\partial_x direction.
%The input image is a 2D image denoted by img,
% theta is the orientation angle of the derivatives, 

% Inputs
%  img: input image
%  theta: orientation angle of left-invariant derivatives
%  npx,npy: number of points in x and y directions respectively
%  dx,dy: spatial steps in the x,y directions.
%  n: order of derivative. E.g., for d^2/dx^2 n = 2.
%  ooa: specify the order of accuracy of the FD scheme. ooa should be an
%      even number to have central difference.
%
%  Outputs:
%  dXi, dEta: (npx*npy)*(npx*npy) sparse matrices

%% Compute the left-invariant derivatives of the 3D lifted image
Df = FDMatrixF(ntheta,npf,df,n,ooa); 
dF = zeros(npx,npy,ntheta, npf); % initialize frequency array
liftedImg = permute(liftedImg, [3 4 1 2]);
for i=1:npx % for each row and column, compute the derivative in the frequency dimension
    for j=1:npy
    oneXindex = liftedImg(:,:,i,j); 
    dfundfNum = Df*oneXindex(:);
    dF(i,j,:,:) = reshape(dfundfNum,ntheta,npf);
    end
end


end