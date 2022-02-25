function kernelFinal = interpolationKernel(vector, intOrder)
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


center = intOrder + 1;
emptyKernel = zeros(2*center-1, 2*center-1);

x = transpose(1:2*center-1);
y = transpose(1:2*center-1);
[X,Y] = meshgrid(x,y);
x = transpose(X(1:end));
y = transpose(Y(1:end));

kernelFinal = emptyKernel;


for i = 1:2*center-1 
    for j = 1:2*center-1   
        kernelOneElement = emptyKernel;
        kernelOneElement(i,j) = 1;
        kernelOneElementInt = scatteredInterpolant(x,y,transpose(kernelOneElement(1:end)));
        kernelFinal(i,j) = kernelOneElementInt(center+vector(1),center+vector(2));  % this is a correlation kernel! To be flipped!         
    end
end


end
