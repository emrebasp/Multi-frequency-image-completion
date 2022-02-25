function kernel = DXiKernelCentral(theta, intOrder)

theta = -theta;
eXi = [cos(theta), sin(theta)];

kernelPlus  = interpolationKernel(eXi, intOrder);
kernelMinus = interpolationKernel(-eXi, intOrder);
kernel      = 0.5* (kernelPlus - kernelMinus);

end