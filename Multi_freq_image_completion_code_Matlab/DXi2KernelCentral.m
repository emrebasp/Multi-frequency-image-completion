function kernel = DXi2KernelCentral(theta, intOrder)

theta = -theta;
eXi = [cos(theta), sin(theta)];

kernelPlus  = interpolationKernel(eXi, intOrder);
kernelCenter = interpolationKernel([0, 0], intOrder);
kernelMinus = interpolationKernel(-eXi, intOrder);
kernel      = kernelPlus - 2*kernelCenter + kernelMinus;

end