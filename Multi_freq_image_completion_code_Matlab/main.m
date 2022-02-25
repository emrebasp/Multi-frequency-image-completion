%%%% Code of the completion algorithm presented in "Multi-frequency image completion via a
%%%% cortically-inspired sub-Riemannian model with frequency and phase", by
%%%% Emre Baspinar, November 2021

%%%% NOTE: This code is based on a single example test image. The parameters
%%%% might need to be changed once a different test image is used.


%% Initialize
clear all;
close all;
clc;


% addpath(strcat(pwd,'/FDMPackage'));

%% Parameters of the Gabor function
sigma     = 2;                                 % : scale of the Gaussian window
minf      = 1.5;                               % : minimum frequency
maxf      = 8;                                 % : maximum frequency
nOfFreq   = 12;                                % : number of frequencies
f         = linspace(minf, maxf, nOfFreq);     % : frequency list of the Gabor filter bank
lambda    = 2*pi*sigma./f;                     % : wavelength, i.e., f=1/lambda
nOftheta  = 32;                                % : number of orientations
theta     = linspace(0,...
2*pi-(2*pi)/nOftheta, nOftheta);               % : orientation angle list for the Gabor function
gamma     = 1;                                 % : spatial aspect ratio (of the x- and y-axis of the Gaussian window)
periodicity = 2*pi;                            % : periodicity of the orientation axis
nOfPhase = 5;                                  % : number of phases  
phase = linspace(0,...                         % : phase list for the Gabor function 
pi/2, nOfPhase);
bandwidth = 0;                                 % : fixed to zero, no control is desired on the spatial frequency bandwidth via this parameter
                                                  
                                            
%% Parameters of the sub-Riemannian diffusion
stepOrientation = (2*pi)/nOftheta;             % : step size in orientation dimension  
stepFreq = (maxf - minf)/(nOfFreq-1);          % : step size in frequency dimension
stepPhi = (pi/2)/(nOfPhase-1);                 % : step size in phase dimension 
imSize  = 256;                                 % : test image size preferably an even number
betaOrientation = nOftheta/sqrt(2*imSize^2);   % : \beta_2 coefficient for unit coherency 
betaPhase = nOfPhase/sqrt(2*imSize^2);         % : \beta_3 coefficient for unit coherency    
betaFreq = nOfFreq/sqrt(2*imSize^2);           % : \beta_4 coefficient for unit coherency
intOrder = 1;                                  % : interpolation order for the B-splines
fTime = 100;                                   % : final time of the diffusion
tStep = 0.1;                                   % : time step of the diffusion

%% Import a test image
% Construct the test image
img = mat2gray(imresize(imread('eyeTest.png'),[imSize, imSize]));
img = normalizer(img(:,:,1));

% Show the test image
figure,
imshow(normalizer(img(8:end-8, 8:end-8)),'InitialMagnification',800);

hold off;

%% Image lifting
% Lift the 2-dim image to the 5-dim sub-Riemannian geometry
for i = 1:nOftheta
    for j = 1:nOfFreq
        for k=1:nOfPhase
     liftRe(:,:,i,j,k) = GaborFilterRe(img, lambda(j), sigma, theta(i), phase(k), gamma, bandwidth); % : real part of the output responses
     liftIm(:,:,i,j,k) = GaborFilterIm(img, lambda(j), sigma, theta(i), phase(k), gamma, bandwidth); % : imaginary part of the output responses
        end
    end
end


%% Create the occluding vertical and horizontal bars at the cortical level
% Configuration for vertical occluding bars
gapX = 2;

xMax0 = 20;
xMin0 = xMax0 -gapX;

xMax1 = 32;
xMin1 = xMax1 -gapX;

xMax2 = 44;
xMin2 = xMax2 -gapX;

xMax3 = 56;
xMin3 = xMax3 -gapX;

xMax4 = 68;
xMin4 = xMax4 -gapX;

xMax5 = 80;
xMin5 = xMax5 -gapX;

xMax6 = 92;
xMin6 = xMax6 -gapX;

xMax7 = 104;
xMin7 = xMax7 -gapX;

xMax8 = 116;
xMin8 = xMax8 -gapX;

xMax9 = 128;
xMin9 = xMax9 -gapX;

xMax10 = 140;
xMin10 = xMax10 -gapX;

xMax11 = 152;
xMin11 = xMax11 -gapX;

xMax12 = 164;
xMin12 = xMax12 -gapX;

xMax13 = 176;
xMin13 = xMax13 -gapX;

xMax14 = 188;
xMin14 = xMax14 -gapX;

xMax15 = 200;
xMin15 = xMax15 -gapX;

xMax16 = 212;
xMin16 = xMax16 -gapX;

xMax17 = 224;
xMin17 = xMax17 -gapX;

xMax18 = 236;
xMin18 = xMax18 -gapX;


% Configuration for horizontal occluding bars
gapY = 2;

yMax0 = 20;
yMin0 = yMax0 - gapY;

yMax1 = 32;
yMin1 = yMax1 - gapY;

yMax2 = 44;
yMin2 = yMax2 - gapY;

yMax3 = 56;
yMin3 = yMax3 - gapY;

yMax4 = 68;
yMin4 = yMax4 - gapY;

yMax5 = 80;
yMin5 = yMax5 - gapY;

yMax6 = 92;
yMin6 = yMax6 - gapY;

yMax7 = 104;
yMin7 = yMax7 - gapY;

yMax8 = 116;
yMin8 = yMax8 - gapY;

yMax9 = 128;
yMin9 = yMax9 - gapY;

yMax10 = 140;
yMin10 = yMax10 - gapY;

yMax11 = 152;
yMin11 = yMax11 - gapY;

yMax12 = 164;
yMin12 = yMax12 - gapY;

yMax13 = 176;
yMin13 = yMax13 - gapY;

yMax14 = 188;
yMin14 = yMax14 - gapY;

yMax15 = 200;
yMin15 = yMax15 - gapY;

yMax16 = 212;
yMin16 = yMax16 - gapY;

yMax17 = 224;
yMin17 = yMax17 - gapY;

yMax18 = 236;
yMin18 = yMax18 - gapY;

% Occluding bars at the cortical level
for i = 1:imSize
    for j = 1:imSize   
        if ((j>=xMin0 &&  j<=xMax0) || (j>=xMin1 &&  j<=xMax1) || (j>=xMin2 &&  j<=xMax2) || (j>=xMin3 &&  j<=xMax3) || (j>=xMin4 &&  j<=xMax4)...
                || (j>=xMin5 &&  j<=xMax5) || (j>=xMin6 &&  j<=xMax6) || (j>=xMin7 &&  j<=xMax7) || (j>=xMin8 &&  j<=xMax8)...
                || (j>=xMin9 &&  j<=xMax9) || (j>=xMin10 &&  j<=xMax10) || (j>=xMin11 &&  j<=xMax11) || (j>=xMin12 &&  j<=xMax12)...
                || (j>=xMin13 &&  j<=xMax13) || (j>=xMin14 &&  j<=xMax14) || (j>=xMin15 &&  j<=xMax15) || (j>=xMin16 &&  j<=xMax16)...
                || (j>=xMin17 &&  j<=xMax17) || (j>=xMin18 &&  j<=xMax18) || (i>=yMin0 && i<=yMax0)...
                || (i>=yMin1 && i<=yMax1) || (i>=yMin2 && i<=yMax2) || (i>=yMin3 && i<=yMax3) || (i>=yMin4 && i<=yMax4)...
                || (i>=yMin5 && i<=yMax5) || (i>=yMin6 && i<=yMax6) || (i>=yMin7 && i<=yMax7) || (i>=yMin8 && i<=yMax8)...
                || (i>=yMin9 && i<=yMax9) || (i>=yMin10 && i<=yMax10) || (i>=yMin11 && i<=yMax11) || (i>=yMin12 && i<=yMax12)...
                || (i>=yMin13 && i<=yMax13) || (i>=yMin14 && i<=yMax14) || (i>=yMin15 && i<=yMax15) || (i>=yMin16 && i<=yMax16)...
                || (i>=yMin17 && i<=yMax17) || (i>=yMin18 && i<=yMax18)) 
            
            liftRe(i,j,:,:,:) = 0; 
            liftIm(i,j,:,:,:) = 0;
        end
    end
end

% Blur a bit the lifted image with occluding bars
D = 5; sigma = 0.25; fsize =  5 * sigma;
gaussFilt = double(NDgaussian(D,fsize,sigma));
liftRe = imfilter(liftRe, gaussFilt);
liftIm = imfilter(liftIm, gaussFilt);
tstImg = liftRe + 1i*liftIm;

% Show the projection of the lifting with occluding bars
liftedProjected = normalizer(sum(liftRe, [3 4 5]));
figure,
imshow(liftedProjected(8:end-8, 8:end-8),'InitialMagnification',800);
hold off;

%% Sub-Riemannian diffusion for completion in the model geometry

allIncrementsInTime = [];
for i=0:tStep:fTime
    
    % Compute all the increments for the corresponding time instant
    incrementXi          = zeros(size(tstImg));
    incrementEta         = zeros(size(tstImg));
    incrementOrientation = zeros(size(tstImg));
    incrementFreq = zeros(size(tstImg));
    incrementPhiPhi = DPhi2Central(tstImg, stepPhi);
    incrementPhi = DPhiCentral(tstImg, stepPhi);     
    
    for thetaLayer = 1:nOftheta
        for fLayer = 1:nOfFreq 
            for phiLayer = 1:nOfPhase
                incrementXi(:,:,thetaLayer,fLayer,phiLayer)   = conv2(tstImg(:,:,thetaLayer,fLayer,phiLayer), DXi2KernelCentral((thetaLayer-1)*(2*pi)/nOftheta, intOrder),'same');
                incrementEta(:,:,thetaLayer,fLayer, phiLayer)  = conv2(tstImg(:,:,thetaLayer,fLayer,phiLayer), DEta2KernelCentral((thetaLayer-1)*(2*pi)/nOftheta, intOrder),'same');
                incrementEtaPhi(:,:,thetaLayer,fLayer, phiLayer)  = 2*f(fLayer)*conv2(incrementPhi(:,:,thetaLayer,fLayer,phiLayer), DEta2KernelCentral(theta(thetaLayer), intOrder),'same'); 
                incrementPhiPhi(:,:,thetaLayer,fLayer, phiLayer) =  (f(fLayer)^2)*incrementPhiPhi(:,:,thetaLayer,fLayer,phiLayer);
            end
        end
    end
    
   incrementOrientation = DTheta2Central(tstImg, stepOrientation);
   incrementFreq        = DFreq2Central(tstImg, stepFreq);
   incrementInTime      = tStep*(incrementXi+betaPhase^2*incrementEta+betaOrientation^2*incrementOrientation+betaFreq^2*incrementFreq...
       +betaPhase^2*incrementEtaPhi + betaPhase^2*incrementPhiPhi);
   
   % Update the lifted image (i.e., the output responses)
   tstImg               = tstImg + incrementInTime; 

   % The boundary conditions are fixed, so restore the source, i.e., the parts which were not occluded at the beginning!     
    for ii = 1:imSize
        for jj = 1:imSize   
            if ~((jj>=xMin0 &&  jj<=xMax0) || (jj>=xMin1 &&  jj<=xMax1) || (jj>=xMin2 &&  jj<=xMax2) || (jj>=xMin3 &&  jj<=xMax3) || (jj>=xMin4 &&  jj<=xMax4)...
                || (jj>=xMin5 &&  jj<=xMax5) || (jj>=xMin6 &&  jj<=xMax6) || (jj>=xMin7 &&  jj<=xMax7) || (jj>=xMin8 &&  jj<=xMax8)...
                || (jj>=xMin9 &&  jj<=xMax9) || (jj>=xMin10 &&  jj<=xMax10) || (jj>=xMin11 &&  jj<=xMax11) || (jj>=xMin12 &&  jj<=xMax12)...
                || (jj>=xMin13 &&  jj<=xMax13) || (jj>=xMin14 &&  jj<=xMax14) || (jj>=xMin15 &&  jj<=xMax15)...
                || (jj>=xMin16 &&  jj<=xMax16) || (jj>=xMin17 &&  jj<=xMax17) || (jj>=xMin18 &&  jj<=xMax18)...
                || (ii>=yMin0 && ii<=yMax0) || (ii>=yMin1 && ii<=yMax1) || (ii>=yMin2 && ii<=yMax2) || (ii>=yMin3 && ii<=yMax3) || (ii>=yMin4 && ii<=yMax4)...
                || (ii>=yMin5 && ii<=yMax5) || (ii>=yMin6 && ii<=yMax6) || (ii>=yMin7 && ii<=yMax7) || (ii>=yMin8 && ii<=yMax8)...
                || (ii>=yMin9 && ii<=yMax9) || (ii>=yMin10 && ii<=yMax10) || (ii>=yMin11 && ii<=yMax11) || (ii>=yMin12 && ii<=yMax12)...
                || (ii>=yMin10 && ii<=yMax10) || (ii>=yMin11 && ii<=yMax11) || (ii>=yMin12 && ii<=yMax12) || (ii>=yMin13 && ii<=yMax13)...
                || (ii>=yMin14 && ii<=yMax14) || (ii>=yMin15 && ii<=yMax15) || (ii>=yMin16 && ii<=yMax16) || (ii>=yMin17 && ii<=yMax17)...
                || (ii>=yMin18 && ii<=yMax18)) 
            
                tstImg(ii,jj,:,:,:) = liftRe(ii,jj,:,:,:) + 1i*liftIm(ii,jj,:,:,:); 
            end
        end
    end
    
    i % keep the track of the time 
    
end


%% Inverse Gabor to come back to the 2-dim image plane
for i = 1:nOftheta
    for j = 1:nOfFreq
        for k = 1:nOfPhase
            invLift(:,:,i,j,k) = sqrt(f(j))*(GaborFilterRe(tstImg(:,:,i,j,k), lambda(j), sigma, theta(i), -phase(k), gamma, bandwidth)+...
            1i*GaborFilterIm(tstImg(:,:,i,j,k), lambda(j), sigma, theta(i), -phase(k), gamma, bandwidth));
        end
    end
end

result = sum(invLift,[3 4 5]);
result = result2(8:end-8, 8:end-8);

% Show the completed image
figure,
imshow(normalizer(result), [0 0.8] ,'InitialMagnification',800);
% title('test image with spike')
hold off;




