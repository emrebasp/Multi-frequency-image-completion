function Dphiphi = DPhi2Central(liftedArray, stepPhi)

RotateForward  = cat(5,liftedArray(:,:,:,:,end), liftedArray(:,:,:,:,1:end-1));
RotateBackward = cat(5,liftedArray(:,:,:,:,2:end), liftedArray(:,:,:,:,1));
Central = liftedArray;
Dphiphi = (1/stepPhi^2) * (RotateForward - 2*Central + RotateBackward); 


end
