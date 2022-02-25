function Dff = DFreq2Central(liftedArray, stepFreq)

RotateForward  = cat(4,liftedArray(:,:,:,end,:), liftedArray(:,:,:,1:end-1,:));
RotateBackward = cat(4,liftedArray(:,:,:,2:end,:), liftedArray(:,:,:,1,:));
Central = liftedArray;
Dff = (1/stepFreq^2) * (RotateForward - 2*Central + RotateBackward); 
 


end
