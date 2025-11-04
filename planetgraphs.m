% orbit params
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

% initial state
x0 = 8; y0 = 1; 
dxdt0 = 0.1; dydt0 = 0;
V0 = [x0, y0, dxdt0, dydt0];

% anonymize func to pass as input later
my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);

BT = get_BT("Fehlberg");

p = 6;

tspan = [1, 3]; h_ref = 0.01;
err_des = 0.01;

[t_list, X_list, h_avg, total_evals, step_fail_rate] = explicit_RK_variable_step_integration(my_rate, tspan, V0, h_ref, BT, p, err_des);

Position = X_list(1:2,:);
Velocity = X_list(3:4,:);

h_list = diff(t_list);

x = Position(1,:);
y = Position(2,:);
r = sqrt( x.^2 + y.^2);
%% Plotting

%For a single desired error threshold, plot the position and velocity as a function of time (position vs. time and velocity vs. time should be on separate axes). Do this with both lines and dots
figure(1)
subplot(1,2,1)
plot(t_list, Position,".-", MarkerSize=15)
title("Planetary Position vs Time")
xlabel("Time"); ylabel("Position")
legend("X Pos","Y Pos")

subplot(1,2,2)
plot(t_list, Velocity,".-", MarkerSize=15)
title("Planetary Velocity vs Time")
xlabel("Time"); ylabel("Velocity")
legend("X Vel","Y Vel")

%For a single desired error threshold, generate a scatter plot of the step size as a function of the distance between the planet and the sun, r = px2 + y2. 
figure(2)
plot(r(:,2:end),h_list,".", MarkerSize=15)
title("Step Size vs Distance Between Planet and Sun (r)")
ylabel("Step Size"); xlabel("Distance")
ylim([-0.2,2])

%% Butcher Tableaus
function BT_struct = get_BT(method_name)

    % initialize structure
    BT_struct = struct();

    % iterate through method names and construct corresponding BT
    if method_name == "Dormand Prince"
        
        BT_struct.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
        BT_struct.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0; ...
        5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
        BT_struct.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0; ...
        3/40, 9/40, 0, 0, 0, 0,0; ...
        44/45, -56/15, 32/9, 0, 0, 0,0; ...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0; ...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0; ...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    elseif method_name == "Fehlberg"

        BT_struct.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
        BT_struct.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55; ...
        25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
        BT_struct.A = [0,0,0,0,0,0; ...
        1/4, 0,0,0,0,0; ...
        3/32, 9/32, 0,0,0,0; ...
        1932/2197, -7200/2197, 7296/2197, 0,0,0; ...
        439/216, -8, 3680/513, -845/4104, 0,0; ...
        -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];

    elseif method_name == "Heun Euler"

        BT_struct.C = [0,1];
        BT_struct.B = [1/2,1/2;1,0];
        BT_struct.A = [0,0;1,0];

    elseif method_name == "Fehlberg RK1"

        BT_struct.C = [0,1/2,1];
        BT_struct.B = [1/512, 255/256, 1/512; ...
        1/256, 255/256, 0];
        BT_struct.A = [0,0,0;1/2,0,0;1/256,255/256,0];
    
    elseif method_name == "Bogacki"

        BT_struct.C = [0,1/2, 3/4, 1];
        BT_struct.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
        BT_struct.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];

    else
        disp("Not a valid method name")
        return 
    end

end
