function [] = init__plot(tit, t, v_true,V0, alpha, theta, q)
%plotter for dutch roll
%just type dr_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%tit = titel, t = time,v_true = true air speed, V0, alpha = angle of attack, theta, pitch, q = pitch rate, 


figure()
sgtitle(tit)

subplot(4,2,1);
plot(t, v_true + V0)
title('True Airspeed V_{TAS} versus Time t')
ylabel('V_{TAS} [m/s]')
xlabel('t [s]')

subplot(4,2,3);
plot(t,alpha)
title('Angle of Attack \alpha versus Time t')
ylabel('\alpha [rad]')
xlabel('t [s]')

subplot(4,2,5);
plot(t,theta)
title('Pitch Angle \theta versus Time t')
ylabel('\theta [rad]')
xlabel('t [s]')

subplot(4,2,7);
plot(t,q)
title('Pitch Rate q versus Time t')
ylabel('q [rad/s]')
xlabel('t [s]')

subplot(4,2,[2,4])
plot(t,alpha)
ylabel('\alpha [rad]')
xlabel('t [s]')
xlim([0,6])

subplot(4,2,[6,8])
plot(t,q)
ylabel('q [rad/s]')
xlabel('t [s]')
xlim([0,6])

end