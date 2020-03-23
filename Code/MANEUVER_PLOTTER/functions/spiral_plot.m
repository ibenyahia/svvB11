function [] = spiral_plot(p_dat, roll_dat, r_dat,delta_a_dat, delta_r_dat, t_dat, p_mod, roll_mod, r_mod, t_mod)
%plotter for spiral
%just type spiral_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%p = roll rate, roll = roll, r = yaw rate, delta_e = elevator def., delta_r =ruder def., t = time


figure;
sgtitle('Spiral (Flight Data)')

subplot(4,1,1);
plot(t_dat,roll_dat,'r',t_mod,roll_mod,'b')
title('Roll angle \phi versus Time t')
ylabel('\phi [rad]')
xlabel('t [s]')
ylim([-1,3])
xlim([0,100])
legend('Flight Data','Model Data')

subplot(4,1,2);
plot(t_dat,p_dat,'r',t_mod,p_mod,'b')
title('Roll rate p versus Time t')
ylabel('p [rad/s]')
xlabel('t [s]')
xlim([0,100])
ylim([-0.1,0.35])


subplot(4,1,3);
plot(t_dat,r_dat,'r',t_mod,r_mod,'b')
title('Yaw rate r versus Time t')
ylabel('r [rad/s]')
xlabel('t [s]')
xlim([0,100])
ylim([-0.1,0.4])

subplot(4,1,4);
hold on
title('Rudder (\delta_r) and Aileron (\delta_a) Deflection versus Time t')
plot(t_dat,delta_r_dat,'g')
plot(t_dat,delta_a_dat,'k')
ylabel('[rad]')
xlabel('t [s]')
xlim([0,100])
ylim([-0.01,0.01])
legend('\delta_r','\delta_a')
hold off


end