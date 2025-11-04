% Run the fixed step integrator using the same Runge-Kutta method that you 
% used for your adaptive step size method across a range of desired step sizes. 

% For your fixed step integrator, remember to use the same row of B in the 
% Butcher tableau that you use in the adaptive step integrator to compute Xn+1.

% Make sure to record the global truncation error, the average step size, and 
% the number of function evaluations.

%% Set up system

clear; clc;

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

%%

p = 5; n = 20;

tspan = [1, 3]; h_ref = 0.05;
h_des = logspace(-6, 2, n);
err_des = logspace(-2, 1, n);

% [h_avg, fcn_evals]
var_nums = zeros(n, 2);
var_g_err = zeros(4, n);
fixed_nums = zeros(n, 2);
fixed_g_err = zeros(4, n);



% fixed step 
for i = 1:n

    [~, X_list, h_avg, total_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_des(i), BT);

    X_final_exact = compute_planetary_motion(tspan(2), V0, orbit_params);
    g_error = abs(X_final_exact - X_list(:, end));

    % store output
    fixed_nums(i, :) = [h_avg, total_evals];
    fixed_g_err(:, i) = g_error;
end



% variable step
for i = 1:n

    [~, X_list, h_avg, total_evals, ~] = explicit_RK_variable_step_integration ...
    (my_rate, tspan, V0, h_ref, BT, p, err_des(i));

    X_final_exact = compute_planetary_motion(tspan(2), V0, orbit_params);
    g_error = abs(X_final_exact - X_list(:, end));

    % store output
    var_nums(i, :) = [h_avg, total_evals];
    var_g_err(:, i) = g_error;
end




%% 

% plot results
figure(1)
% g_err vs avg step size
loglog(var_nums(:,1), var_g_err(1,:), ".-", MarkerSize=15, DisplayName="Adaptive Step Size")
% hold on;
% loglog(fixed_nums(:,1), fixed_g_err(1,:), ".-", DisplayName="Fixed Step Size")

% label plots
title("Global Truncation Error vs Step Size")
xlabel("Average Step Size"); ylabel("Error")
legend()

hold off; figure(2)
loglog(var_nums(:,2), var_g_err(1,:), ".-", MarkerSize=15, DisplayName="Adaptive Step Size"); 
% hold on;
% loglog(fixed_nums(:,2), fixed_g_err(1,:), ".-", DisplayName="Fixed Step Size")

% label plots
title("Global Truncation Error vs Step Size")
xlabel("Num. Function Evaluations"); ylabel("Error")
legend()


%%

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
