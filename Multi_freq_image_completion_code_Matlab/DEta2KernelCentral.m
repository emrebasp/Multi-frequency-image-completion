function kernel = DEta2KernelCentral(theta, intOrder)

theta = pi-theta;
eEta  = [-sin(theta), cos(theta)];

kernelPlus  = interpolationKernel(eEta, intOrder);
kernelCenter = interpolationKernel([0, 0], intOrder);
kernelMinus = interpolationKernel(-eEta, intOrder);
kernel      = kernelPlus - 2*kernelCenter + kernelMinus;

end