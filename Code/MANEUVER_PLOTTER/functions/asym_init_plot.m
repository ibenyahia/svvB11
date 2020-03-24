function [] = asym_init_plot(tit, t, beta, phi, p, r)
%plotter for dutch roll
%just type dr_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%tit = titel, t = time, beta = side slip, phi = roll, p = roll rate, r = yaw rate, 


figure()
sgtitle(tit)

subplot(2,2,1);
plot(t, beta)
title('Side slip \beta versus Time t')
ylabel('\beta [rad]')
xlabel('t [s]')

subplot(2,2,2);
plot(t,phi)
title('Roll angle \phi versus Time t')
ylabel('\phi [rad]')
xlabel('t [s]')

subplot(2,2,3);
plot(t,r)
title('Yaw rate r versus Time t')
ylabel('r [rad/s]')
xlabel('t [s]')

subplot(2,2,4);
plot(t,p)
title('Roll rate p Rate versus Time t')
ylabel('p [rad/s]')
xlabel('t [s]')

end