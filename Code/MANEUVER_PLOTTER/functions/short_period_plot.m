function [] = short_period_plot(alpha_dat,theta_dat, delta_e_dat, q_dat, t_dat, alpha_mod,q_mod, t_mod,theta_mod)
%plotter for short period eigenmotion
%just type short_period_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%theta = pitch ,delta = elevator pitch, q = pitch rate , t = time 
figure;
sgtitle('Short Period Improved')

subplot(2,2,2);
plot(t_dat,theta_dat,'r',t_mod,theta_mod,'b')
title(' Pitch \theta versus Time t')
ylabel('\theta [rad]')
xlabel('t [s]')
legend('Flight Data','Model Data')

subplot(2,2,4);
plot(t_dat,q_dat,'r',t_mod,q_mod,'b')
title('Pitch rate q versus Time t')
ylabel('q [rad/s]')
xlabel('t [s]')

subplot(2,2,1);
plot(t_dat, delta_e_dat,'g')
title('Elevator deflection \delta_e versus Time t')
ylabel('\delta_e [rad]')
xlabel('t [s]')

subplot(2,2,3);
plot(t_dat,alpha_dat,'r',t_mod,alpha_mod,'b')
title('Angle of attack \alpha versus Time t')
ylabel('\alpha [rad]')
xlabel('t [s]')

end