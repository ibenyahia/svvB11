function [] = short_period_plot(theta_dat, delta_dat,theta_dt_dat, t_dat, theta_mod, delta_mod,theta_dt_mod, t_mod)
%plotter for short period eigenmotion
%just type short_period_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%theta = pitch ,delta = elevator pitch, theta_dt = pitch rate , t = time 
figure;
sgtitle('Short Period')

subplot(2,2,1);
plot(t_dat,theta_dat,'r',t_mod,theta_mod,'b')
title('pitch vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,3);
plot(t_dat,theta_dt_dat,'r',t_mod,theta_dt_mod,'b')
title('pitch rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,2);
plot(t_dat,delta_dat,'r',t_mod,delta_mod,'b')
title('delta_e vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')


end