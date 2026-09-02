function [rv_FRT, rv0_FRT_SATs, rv0_NRHO_south_SATs, rv0_L4_vertical_SATs, rv0_L5_vertical_SATs, tf_FRT] = ...
    generate_initial_rv_data(number_of_timesteps, number_of_FRT_SATs, number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs)

global mu 
% Free Return Trajectory
rv0_FRT = [-0.024087817019224;
           -0.013199742313013;
            0;
            7.725332506766247;
           -7.042775379786875;
            0];

tf_FRT = 1.304372141714530;
t_FRT_vector = (0:tf_FRT/number_of_timesteps:tf_FRT)';

odeopt = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_FRT, rv_FRT] = ode45(@CRTBP, t_FRT_vector, rv0_FRT, odeopt);


%% NRHO-South #742

rv0_NRHO_south = [1.0118016478467127E+0
                  2.7862723172731632E-27
                  -1.7388786181524768E-1	
                  -6.5977076420792335E-13	
                  -7.9865444979787661E-2	
                  -9.5808807140677979E-12];

NRHO_south_period = 1.3739476422245205E+0;
t_NRHO_vector = 0:NRHO_south_period/number_of_timesteps:NRHO_south_period;


[t_NRHO_south, rv_NRHO_south] = ode45(@CRTBP, t_NRHO_vector, rv0_NRHO_south, odeopt);


% L4 / L5 Vertical Orbits

%% L4 Vertical #1020

rv0_L4_vertical = [6.3891964038363835E-1	
                   4.2544999999999999E-1	
                   6.3338764777407142E-1	
                   -3.1073088790628428E-1	
                   -4.8651604983496066E-1	
                   6.4806934253088289E-1];

L4_period = 6.2976761954249154E+0;
t_L4_vector = 0:L4_period/number_of_timesteps:L4_period;
[t_L4, rv_L4] = ode45(@CRTBP, t_L4_vector,rv0_L4_vertical, odeopt);



%% L5 Vertical #1020

rv0_L5_vertical = [6.3891964038363835E-1	
                   -4.2544999999999999E-1	
                   6.3338764777407142E-1	
                   3.1073088790628428E-1	
                   -4.8651604983496066E-1	
                   -6.4806934253088289E-1];



[t_L5, rv_L5] = ode45(@CRTBP, t_L4_vector,rv0_L5_vertical, odeopt);


rv0_FRT_SATs = rv_FRT(1:number_of_timesteps/number_of_FRT_SATs:(number_of_timesteps -number_of_timesteps/number_of_FRT_SATs+1) ,1:6);
rv0_NRHO_south_SATs  = rv_NRHO_south(1:number_of_timesteps/number_of_NRHO_SATs:(number_of_timesteps -number_of_timesteps/number_of_NRHO_SATs+1) ,1:6);
rv0_L4_vertical_SATs = rv_L4(1:number_of_timesteps/number_of_L4_SATs:(number_of_timesteps -number_of_timesteps/number_of_L4_SATs+1) ,1:6);
rv0_L5_vertical_SATs = rv_L5(1:number_of_timesteps/number_of_L5_SATs:(number_of_timesteps -number_of_timesteps/number_of_L5_SATs+1) ,1:6);

figure;
scatter3(rv0_FRT_SATs(:,1),rv0_FRT_SATs(:,2),rv0_FRT_SATs(:,3),"*",'r')
hold on
scatter3(rv0_NRHO_south_SATs(:,1),rv0_NRHO_south_SATs(:,2),rv0_NRHO_south_SATs(:,3),'*','m')
scatter3(rv0_L4_vertical_SATs(:,1),rv0_L4_vertical_SATs(:,2),rv0_L4_vertical_SATs(:,3),'*','g')
scatter3(rv0_L5_vertical_SATs(:,1),rv0_L5_vertical_SATs(:,2),rv0_L5_vertical_SATs(:,3),'*','b')
scatter3(-mu,0,0,'o','filled')
scatter3(1-mu,0,0,'o','filled')

plot3(rv_FRT(:,1),rv_FRT(:,2), rv_FRT(:,3),'r')
plot3(rv_NRHO_south(:,1),rv_NRHO_south(:,2), rv_NRHO_south(:,3),'m')
plot3(rv_L4(:,1),rv_L4(:,2),rv_L4(:,3),'g')
plot3(rv_L5(:,1),rv_L5(:,2),rv_L5(:,3),'b')
hold off

title('FRT satellites and Cislunar PNT constellation','FontSize',15,'FontWeight','bold')
xlabel('x [LU]')
ylabel('y [LU]')
zlabel('z [LU]')
legend('Free Return Trajectory', 'L2 NRHO South', 'L4 Vertical', 'L5 Vertical','Earth','Moon')





end