function [] = spiral_plot(p_dat, roll_dat, r_dat,delta_e_dat, delta_r_dat, t_dat, p_mod, roll_mod, r_mod, t_mod, damper)
%plotter for spiral
%just type spiral_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%p = roll rate, roll = roll, r = yaw rate, delta_e = elevator def., delta_r =ruder def., t = time


figure;
sgtitle('Spiral')

subplot(2,3,2);
plot(t_dat,roll_dat,'r',t_mod,roll_mod,'b')
title('roll vs t')
ylabel('[rad]')
xlabel('[s]')
legend('Flight Data','Model Data')

subplot(2,3,5);
plot(t_dat,p_dat,'r',t_mod,p_mod,'b')
title('roll rate vs t')
ylabel('[rad/s]')
xlabel('[s]')



subplot(2,3,3);
plot(t_dat,r_dat,'r',t_mod,r_mod,'b')
title('yaw rate vs t')
ylabel('[rad/s]')
xlabel('[s]')


subplot(2,3,1);
plot(t_dat,delta_r_dat,'g')
title('delta_r vs t')
ylabel('[rad]')
xlabel('[s]')
legend('Control input')

subplot(2,3,4);
plot(t_dat,delta_e_dat,'g')
title('delta_e vs t')
ylabel('[rad]')
xlabel('[s]')

end