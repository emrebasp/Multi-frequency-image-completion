function dEta = DEta(dx,dy,theta, npx, npy)
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


theta = - theta;
%% Compute the left-invariant derivatives of the 2D image
dEta = -sin(theta)*dx + cos(theta)*dy; 
dEta = reshape(dEta,npx,npy);

end