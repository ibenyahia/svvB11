load('cit_par_bigdataourflight.mat') %File with all the constants and variables

%Symmetric eigenmotions

%Short period motion
A_SPM = 2*muc_SPM*KY2*(2*muc_SPM-CZadot);
B_SPM = -2*muc_SPM*KY2*CZa-(2*muc_SPM+CZq)*Cmadot-(2*muc_SPM-CZadot)*Cmq;
C_SPM = CZa*Cmq-(2*muc_SPM+CZq)*Cma;

Lambda_1s = (-B_SPM+sqrt(-(4*A_SPM*C_SPM)+(B_SPM^2)))/(2*A_SPM);
%Lambda_1sscalar = sqrt(real(Lambda_1s)^2+imag(Lambda_1s)^2)

% use V0 for V???
%P_SPM1 = (2*pi/imag(Lambda_1s))*(c/V0(1))
%T_SPM1 = (lnog1/2)/Lambda_1sscalar)*(c/V0(1))
%D_SPM1 = -real(Lambda_1s)/(Lambda_1sscalar)

Lambda_2s = (-B_SPM-sqrt(-(4*A_SPM*C_SPM)+(B_SPM^2)))/(2*A_SPM);
%Lambda_2sscalar = sqrt(real(Lambda_2s)^2+imag(Lambda_2s)^2)

%P_SPM2 = (2*pi/imag(Lambda_2s))*(c/V0(1))
%T_SPM2 = (log(1/2)/Lambda_2sscalar)*(c/V0(1))
%D_SPM2 = -real(Lambda_2s)/(Lambda_2sscalar)


%Phugoid
A_PHU = 2*muc_PHU*(CZa*Cmq-2*muc_PHU*Cma);
B_PHU = 2*muc_PHU*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa);
C_PHU = CZ0_PHU*(Cmu*CZa-CZu*Cma);

Lambda_3s = (-B_PHU+sqrt(-4*A_PHU*C_PHU+B_PHU^2))/(2*A_PHU);
%Lambda_3sscalar = sqrt(real(Lambda_3s)^2+imag(Lambda_3s)^2)

%P_PHU1 = (2*pi/imag(Lambda_3s))*(c/V0(2))
%T_PHU1 = (log(1/2)/Lambda_3sscalar)*(c/V0(2))
%D_PHU1 = -real(Lambda_3s)/(Lambda_3sscalar)

Lambda_4s = (-B_PHU-sqrt(-4*A_PHU*C_PHU+B_PHU^2))/(2*A_PHU);
%Lambda_4sscalar = sqrt(real(Lambda_4s)^2+imag(Lambda_4s)^2)

%P_PHU2 = (2*pi/imag(Lambda_4s))*(c/V0(2))
%T_PHU2 = (log(1/2)/Lambda_4sscalar)*(c/V0(2))
%D_PHU2 = -real(Lambda_4s)/(Lambda_4sscalar)




%Asymmetric eigenmotions

%Aperiodic roll motion

Lambda_1a = Clp/(4*mub_ARM*KX2);
%Lambda_1ascalar = sqrt(real(Lambda_1a)^2+imag(Lambda_1a)^2)

%P_ARM = (2*pi/imag(Lambda_1a))*(c/V0(3))
%T_ARM = (log(1/2)/Lambda_1ascalar)*(c/V0(3))
%D_ARM = -real(Lambda_1a)/(Lambda_1ascalar)

%Dutch roll (only use DR not DRDAMP, bcs no equations for that one)
A_DR = 8*(mub_DR^2)*KZ2 ;
B_DR = -2*mub_DR*(Cnr+2*KZ2*CYb) ;
C_DR = (4*mub_DR*Cnb)+(CYb*Cnr) ;

Lambda_2a = (-B_DR+sqrt(-(4*A_DR*C_DR)+(B_DR^2)))/(2*A_DR);
%Lambda_2ascalar = sqrt(real(Lambda_2a)^2+imag(Lambda_2a)^2)

%P_DRM1 = (2*pi/imag(Lambda_2a))*(c/V0(4))
%T_DRM1 = (log(1/2)/Lambda_2ascalar)*(c/V0(4))
%D_DRM1 = -real(Lambda_2a)/(Lambda_2ascalar)

Lambda_3a = (-B_DR-sqrt(-(4*A_DR*C_DR)+(B_DR^2)))/(2*A_DR);
%Lambda_3ascalar = sqrt(real(Lambda_3a)^2+imag(Lambda_3a)^2)

%P_DRM2 = (2*pi/imag(Lambda_3a))*(c/V0(4))
%T_DRM2 = (log(1/2)/Lambda_3ascalar)*(c/V0(4))
%D_DRM2 = -real(Lambda_3a)/(Lambda_3ascalar)

%Spiral motion

Lambda_4a = (2*CL_SM*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub_SM*Cnb)-Cnp*(CYb*Clr+4*mub_SM*Clb));
%Lambda_4ascalar = sqrt(real(Lambda_4a)^2+imag(Lambda_4a)^2)

%P_SM = (2*pi/imag(Lambda_4a))*(c/V0(6))
%T_SM = (log(1/2)/Lambda_4ascalar)*(c/V0(6))
%D_SM = -real(Lambda_4a)/(Lambda_4ascalar)

