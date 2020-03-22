function [] = phugoid_plot(vtas_dat, theta_dat, delta_e_dat, q_dat, t_dat, vtas_mod, theta_mod, q_mod,  t_mod)
%plotter for short phugoid
%just type phugoid_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%vtas = true air spedd, theta = pitch ,delta = elevator pitch, q = pitch rate, t = time 

figure;
sgtitle('Phugoid')

subplot(2,3,[2,3]);
plot(t_dat,theta_dat,'r',t_mod,theta_mod,'b')
title('pitch vs t')
ylabel('[rad]')
xlabel('[s]')
legend('Flight Data','Model Data')

subplot(3,3,[5,6]);
plot(t_dat,q_dat,'r',t_mod,q_mod,'b')
title('pitch rate vs t')
ylabel('[rad]')
xlabel('[s]')


subplot(3,3,4);
plot(t_dat,delta_e_dat,'g')
title('delta_e vs t')
ylabel('[deg]')
xlabel('[s]')
legend('Control input')

subplot(3,3,4);
plot(t_dat,vtas_dat,'r',t_mod,vtas_mod,'b')
title('Speed vs t')
ylabel('[m/s]')
xlabel('[s]')



end