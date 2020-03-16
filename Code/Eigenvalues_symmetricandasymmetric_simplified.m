load('variables.mat') %File with all the constants and variables

%Symmetric eigenmotions

%Short period motion
A_SPM = 2*muc*KY2*(2*muc-CZadot);
B_SPM = -2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq;
C_SPM = CZa*Cmq-(2*muc+Cma)*Cma;

Lambda_1s = (-B_SPM+sqrt(4*A_SPM*C_SPM-B_SPM^2))/(2*A_SPM)

Lambda_2s = (-B_SPM-sqrt(4*A_SPM*C_SPM-B_SPM^2))/(2*A_SPM)


%Phugoid
A_PHM = 2*muc*(CZa*Cmq-2*muc*Cma);
B_PHM = 2*muc*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa);
C_PHM = CZ0*(Cmu*CZa-CZu*Cma);

Lambda_3s = (-B+sqrt(4*A*C-B^2))/(2*A)

Lambda_4s = (-B-sqrt(4*A*C-B^2))/(2*A)



%Asymmetric eigenmotions

%Aperiodic roll motion

Lambda_1a = Clp/(4*mub*KX2)

%Dutch roll

Lambda_2a = 2*(Cnr+2*KZ2*CYb)+sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*Cyb)^2)/(16*mub*KZ2)

Lambda_3a = 2*(Cnr+2*KZ2*CYb)-sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*Cyb)^2)/(16*mub*KZ2)

%Spiral motion

Lambda_4a = (2*CL*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))


