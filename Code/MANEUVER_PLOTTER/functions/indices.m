function [i_0,i_1] = indices(h_0,min_0,sec_0,duar,t)
% this function will return the index of the closest value of the begining time and end time 
%h_0 = hour time stamp of motion begin
%min_0 = minute time stamp
%sec_0 = second time stamp
%duar = duaration of maneuver in second
%t = list off time intervals of all measurements 
t_0 = h_0 * 60^2 + min_0 * 60 + sec_0 ;
t_1 = t_0 + duar ;

[dif_sp_0,i_0] = min(abs(t-t_0)) ;        %gives index of start point
[dif_sp_1,i_1] = min(abs(t-t_1)) ;       %gives index of end point
end 

