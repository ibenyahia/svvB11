function [] = ap_roll_plot(roll_rate_dat, delta_dat,roll_dat, t_dat, roll_rate_mod, delta_mod, roll_mod, t_mod)
%plotter for short aperiodic roll
%just type ap_roll_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%roll_rate = roll rate, roll = roll ,delta = elevator pitch, t = time 

figure;
sgtitle('Aperiodic Roll')

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
plot(t_dat,delta_dat,'r',t_mod,delta_mod,'b')
title('delta_e vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')


end