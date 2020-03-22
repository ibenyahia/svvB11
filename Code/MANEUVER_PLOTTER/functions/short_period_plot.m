function [] = short_period_plot(alpha_dat,theta_dat, delta_dat, q_dat, t_dat, alpha_mod,q_mod, t_mod)
%plotter for short period eigenmotion
%just type short_period_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%theta = pitch ,delta = elevator pitch, q = pitch rate , t = time 
figure;
sgtitle('Short Period')

subplot(2,2,1);
plot(t_dat,theta_dat,'r',t_mod,theta_mod,'b')
title('pitch vs t')
ylabel('[rad]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,3);
plot(t_dat,q_dat,'r',t_mod,q_mod,'b')
title('pitch rate vs t')
ylabel('[rad/s]')
xlabel('[s]')

subplot(2,2,4);
plot(t_dat, delta_dat,'g')
title('delta_e vs t')
ylabel('[rad]')
xlabel('[s]')

subplot(2,2,2);
plot(t_dat,alpha_dat,'r',t_mod,alpha_mod,'b')
title('AoA vs t')
ylabel('[rad]')
xlabel('[sec]')

end