function [] = spiral_plot(roll_rate_dat, h_dat,roll_dat, t_dat, roll_rate_mod, h_mod, roll_mod, t_mod)
%plotter for spiral
%just type spiral_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%roll_rate = roll rate, h = altitude, roll = roll , t = time 


figure;
sgtitle('Spiral')

subplot(2,2,1);
plot(t_dat,roll_dat,'r',t_mod,roll_mod,'b')
title('roll vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,3);
plot(t_dat,roll_rate_dt_dat,'r',t_mod,roll_rate_mod,'b')
title('roll rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,2);
plot(t_dat,h_dat,'r',t_mod,h_mod,'b')
title('altitude vs t')
ylabel('[km]')
xlabel('[sec]')
legend('Flight Data','Model Data')


end