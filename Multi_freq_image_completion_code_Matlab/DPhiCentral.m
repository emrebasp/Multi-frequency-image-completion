function Dphi = DPhiCentral(liftedArray, stepPhi)

RotateForward  = cat(5,liftedArray(:,:,:,:,end), liftedArray(:,:,:,:,1:end-1));
RotateBackward = cat(5,liftedArray(:,:,:,:,2:end), liftedArray(:,:,:,:,1));
Dphi = 1/(2*stepPhi) * (RotateForward - RotateBackward); 


end
