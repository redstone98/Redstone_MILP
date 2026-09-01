% Consider the planar Earth-Moon RTBP, μ=0.01215
% A) Solve for exact locations of L1-L5 (in typical RTBP frame with origin at barycenter). Express results
% in LU.

clear; 

m1 = 398600.4418; 
m2 = 4902.800076; 

global mu 

mu = m2 / (m1 + m2); 

% ------------------------------------------------------------------------
% L2 
x0 = 1.1; 
x2 = fsolve(@L2, x0); 
y2 = 0; 
% ------------------------------------------------------------------------
% L1  
x0 = 0.85; 
x1 = fsolve(@L1, x0); 
y1 = 0; 
% ------------------------------------------------------------------------
% L3 
x0 = -1.45; 
x3 = fsolve(@L3, x0); 
y3 = 0; 
% ------------------------------------------------------------------------
% L4 
x4 = -mu + 1/2; 
y4 = sqrt(3)/2; 
% ------------------------------------------------------------------------
% L5 
x5 = -mu+1/2; 
y5 = -sqrt(3)/2; 


%% solved ICs 

% % L1 north 
% rv0 = [
%   -0.024087817019224
%   -0.013199742313013
%                    0
%    7.734537264743389
%   -7.032951682588838
%                    0
%        ]; 
% 
%     % propagate halo orbit 
%     toler   = 1e-8;
%     tf = 1.34504565372; 
%     options = odeset('reltol', toler, 'abstol', 1e-8); 
%     [t, rv_FRT] = ode45(@CRTBP, [0 tf], rv0, options); 
 
%% Initial condition
rv0_FRT = [
  -0.024087817019224
  -0.013199742313013
                   0
   7.734537264743389
  -7.032951682588838
                   0
];

tf0 = 1.34504565372;

toler = 1e-10;
odeopt = odeset('RelTol',toler,'AbsTol',1e-10);

%% Target symmetric point
r_target = [rv0_FRT(1); -rv0_FRT(2)];

%% Decision variables
% q = [vx0 vy0 tf]

q0 = [rv0_FRT(4); rv0_FRT(5); tf0];

%% Bounds
dv = 0.5;
dt = 0.2;

lb = q0 + [-dv; -dv; -dt];
ub = q0 + [ dv;  dv;  dt];

%% fmincon options
opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',500, ...
    'OptimalityTolerance',1e-12, ...
    'StepTolerance',1e-12');

%% Solve
q_sol = fmincon( ...
    @(q) symmetry_cost(q, rv0_FRT, r_target, odeopt), ...
    q0, ...
    [],[],[],[], ...
    lb, ub, ...
    [], ...
    opts);

%% Updated state
rv0_new = rv0_FRT;

rv0_new(4) = q_sol(1);
rv0_new(5) = q_sol(2);

tf_new = q_sol(3);

%% Final propagation
[t, rv_FRT] = ode45(@CRTBP, [0;tf_new], rv0_new, odeopt);

r_final = rv_FRT(end,1:2).';

disp('Corrected initial condition:')
disp(rv0_new)

disp('Corrected tf:')
disp(tf_new)

disp('Final position:')
disp(r_final)

disp('Target:')
disp(r_target)

disp('Residual:')
disp(r_final-r_target)
%% NRHO-South #742

rv0_NRHO_south = [1.0118016478467127E+0
                  2.7862723172731632E-27
                  -1.7388786181524768E-1	
                  -6.5977076420792335E-13	
                  -7.9865444979787661E-2	
                  -9.5808807140677979E-12];

NRHO_south_period = 1.3739476422245205E+0;

[t_NRHO_south, rv_NRHO_south] = ode45(@CRTBP, [0 10*NRHO_south_period], rv0_NRHO_south, odeopt);


%% NRHO-North #742

rv0_NRHO_north = [1.0113254829162490E+0	
                  -3.4655306799190984E-28	
                  1.7343215557041181E-1	
                  -6.8885128080822615E-13	
                  -7.8717210400801721E-2	
                  1.0116298133884205E-11];

NRHO_north_period = 1.3672906570047736E+0;

[t_NRHO_north, rv_NRHO_north] = ode45(@CRTBP, [0 10*NRHO_north_period], rv0_NRHO_north, odeopt);

%% L4 Vertical #1020

rv0_L4_vertical = [6.3891964038363835E-1	
                   4.2544999999999999E-1	
                   6.3338764777407142E-1	
                   -3.1073088790628428E-1	
                   -4.8651604983496066E-1	
                   6.4806934253088289E-1];

L4_period = 6.2976761954249154E+0;


[t_L4, rv_L4] = ode45(@CRTBP, [0 10 * L4_period],rv0_L4_vertical, odeopt);



%% L5 Vertical #1020

rv0_L5_vertical = [6.3891964038363835E-1	
                   -4.2544999999999999E-1	
                   6.3338764777407142E-1	
                   3.1073088790628428E-1	
                   -4.8651604983496066E-1	
                   -6.4806934253088289E-1];

L5_period = 6.2976761954249154E+0;


[t_L5, rv_L5] = ode45(@CRTBP, [0 10*L5_period],rv0_L5_vertical, odeopt);


%% L4 Short Period #1050

rv0_L4_short = [4.8784941344943100E-1	
                -1.0404208063235154E+0	
                7.1991551496244772E-25	
                -2.6678436944953493E-1	
                -1.1575102652787685E-1	
                -9.9768091232302432E-26];

L4_short_period = 2.9109498144129208E+1;


[t_L4_short, rv_L4_short] = ode45(@CRTBP, [0 10*L4_short_period],rv0_L4_short, odeopt);




%% L5 Short Period #1050

rv0_L5_short = [4.8784941344943100E-1	
                   1.0404208003037976E+0	
                   -2.4452718823520319E-24	
                   2.6678436047453685E-1	
                   -1.1575102236407500E-1	
                   2.3627799876915897E-25];

L5_short_period = 2.9109498148862315E+1;


[t_L5_short, rv_L5_short] = ode45(@CRTBP, [0 10*L5_short_period],rv0_L5_short, odeopt);

%% Plot
figure;
plot(rv_FRT(:,1), rv_FRT(:,2),'LineWidth',1.5)
hold on

plot(rv0_new(1), rv0_new(2), ...
    'ro','MarkerFaceColor','r')

plot(r_final(1), r_final(2), ...
    'bo','MarkerFaceColor','b')

plot(r_target(1), r_target(2), ...
    'kx','LineWidth',2,'MarkerSize',10)

axis equal
grid on

xlabel('x')
ylabel('y')

legend('Trajectory','Initial','Final','Target')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COST FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = symmetry_cost(q, rv0, r_target, odeopt)

    rv_init = rv0;

    % Adjust vx vy
    rv_init(4) = q(1);
    rv_init(5) = q(2);

    tf = q(3);

    [~, rv] = ode45(@CRTBP, [0 tf], rv_init, odeopt);

    r_final = rv(end,1:2).';

    err = r_final - r_target;

    J = err.'*err;

end

    % plot 
    
% create figure 
name = 'CRTBP: FRT PNT Servicing'; 
figure('name', name); 
    scatter3(x1, y1, 0, 'p','DisplayName', 'L1'); hold on; 
    scatter3(x2, y2, 0, 'p', 'DisplayName', 'L2'); 
    scatter3(x4, y4, 0, 'p', 'DisplayName', 'L4'); 
    scatter3(x5, y5, 0, 'p', 'DisplayName', 'L5'); 
    scatter3(1 - mu, 0, 0, 20, 'filled', 'DisplayName', 'Moon'); 
    scatter3(-mu, 0, 0, 50, 'filled', 'DisplayName', 'Earth'); 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 
    title(name) 
    
    % L1 N 
    plot_orbit_0f(rv_FRT);
    plot_orbit_0f(rv_NRHO_south);
    plot_orbit_0f(rv_NRHO_north);
    plot_orbit_0f(rv_L4);
    plot_orbit_0f(rv_L5);
    plot_orbit_0f(rv_L4_short);
    plot_orbit_0f(rv_L5_short);
    
    legend('L1', 'L2', 'L4', 'L5','Moon', 'Earth'); 
    axis('equal')
%% subfunctions 


