function Dthetatheta = DTheta2Central(liftedArray, stepTheta)

RotateForward  = cat(3,liftedArray(:,:,end,:,:), liftedArray(:,:,1:end-1,:,:));
RotateBackward = cat(3,liftedArray(:,:,2:end,:,:), liftedArray(:,:,1,:,:));
Central = liftedArray;
Dthetatheta = (1/stepTheta^2) * (RotateForward - 2*Central + RotateBackward); 


end
