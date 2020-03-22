function [] = dr_plot(roll_rate_dat, roll_dat, yaw_rate_dat, t_dat, roll_rate_mod, roll_mod, yaw_rate_mod, t_mod, damper)
%plotter for dutch roll
%just type dr_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%roll_rate = roll rate, roll = roll, yaw = yaw, yaw_rate = yaw rate, t = time 


figure;
if damper == 1
    sgtitle('Damped Dutch Roll')
else
    sgtitle('Undamped Dutch Roll')
end

subplot(2,2,1);
plot(t_dat,roll_dat,'r',t_mod,roll_mod,'b')
title('roll vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,3);
plot(t_dat,roll_rate_dat,'r',t_mod,roll_rate_mod,'b')
title('roll rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')
legend('Flight Data','Model Data')


subplot(2,2,2);
plot(t_dat,yaw_rate_dat,'r',t_mod,yaw_rate_mod,'b')
title('yaw rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')
legend('Flight Data','Model Data')

end