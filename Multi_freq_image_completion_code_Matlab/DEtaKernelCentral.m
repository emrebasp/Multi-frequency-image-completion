function kernel = DEtaKernelCentral(theta, intOrder)

theta = pi-theta;

eEta  = [-sin(theta), cos(theta)];

kernelPlus  = interpolationKernel(eEta, intOrder);
kernelMinus = interpolationKernel(-eEta, intOrder);
kernel      = 0.5* (kernelPlus - kernelMinus);

end 