function [Dtheta, Df] = DThetaFreqKernelCentral(liftedArray)

[Dx Dy Dtheta Df] = gradient(liftedArray);
end