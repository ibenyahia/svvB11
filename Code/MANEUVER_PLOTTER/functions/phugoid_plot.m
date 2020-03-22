function [] = phugoid_plot(V_dat, h_dat,theta_dat,delta_dat, t_dat, V_mod, h_mod, theta_mod, delta_mod,  t_mod)
%plotter for short phugoid
%just type phugoid_plot and plug in your data zack beautiful graph
%_dat for flight/ref data, _mod for modeled data
%V = spedd, h = altutude theta = pitch ,delta = elevator pitch, t = time 

figure;
sgtitle('Phugoid')

subplot(2,2,1);
plot(t_dat,theta_dat,'r',t_mod,theta_mod,'b')
title('pitch vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,3);
plot(t_dat,h_dat,'r',t_mod,h_mod,'b')
title('altitude vs t')
ylabel('[km]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,2);
plot(t_dat,delta_dat,'r',t_mod,delta_mod,'b')
title('delta_e vs t')
ylabel('[deg]')
xlabel('[sec]')
legend('Flight Data','Model Data')

subplot(2,2,4);
plot(t_dat,V_dat,'r',t_mod,V_mod,'b')
title('Speed vs t')
ylabel('[m/s]')
xlabel('[sec]')
legend('Flight Data','Model Data')


end