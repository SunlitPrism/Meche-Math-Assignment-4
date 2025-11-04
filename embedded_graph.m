% for selected embedded method

% How do the local truncation errors of XB1 and XB2 scale with h (find the corresponding values of p). 
% Which of the two has better accuracy? How does |XB1 − XB2| scale with h?

% Plot the local truncation errors of XB1 and XB2 as a function of their difference, |XB1 − XB2|. How do the errors scale with |XB1 − XB2|? Is |XB1 − XB2| a good proxy for the error (or at least, does it provide a good upper bound on the error?).
%% Housekeeping
clear; clc; clf;
C = rgb2hex(orderedcolors("gem"));

%% Initialize system

% selected method: Fehlberg
BT = get_BT("Fehlberg");

% orbit params
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

% initial state
x0 = 8; y0 = 1; 
dxdt0 = 0.1; dydt0=0;
V0 = [x0, y0, dxdt0, dydt0];

% anonymize func to pass as input later
my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);


%% Generate a loglog plot of the LTE for XB1 and XB2 as a function of the step size, h.
% On the same axes, plot |XB1 − XB2| as a function of h, as well as the difference |f (tref + h) − f (t)|.

% intialize some variables
n = 100; h_range = logspace(-5, 1, n); 
method_name = "Fehlberg"; t_ref = 0.331; 
XA = compute_planetary_motion(t_ref, V0, orbit_params);

% compute local truncation error
LTE_XB1 = get_local_tr_err(my_rate, XA, V0, t_ref, 1, h_range, BT, orbit_params);
LTE_XB2 = get_local_tr_err(my_rate, XA, V0, t_ref, 2, h_range, BT, orbit_params);
XB_diff = abs(LTE_XB1 - LTE_XB2);

% compute reference difference
ref_list = zeros(length(XA), n);
for i = 1:n
    V_list = compute_planetary_motion([t_ref, t_ref+h_range(i)], V0, orbit_params);
    ref_list(:, i) = abs(diff(V_list));
end

% plot results
figure(1);
loglog(h_range, ref_list(1, :), "-", DisplayName="Ref. list"); hold on;
loglog(h_range, LTE_XB1, ".-", DisplayName=(method_name + ": XB1"));
loglog(h_range, LTE_XB2, ".-", DisplayName=(method_name + ": XB2"));
loglog(h_range, XB_diff, ".", DisplayName=(method_name + ": |XB1-XB2|"));

% label plots
title("Fehlberg Local Truncation Error")
xlabel("Step Size"); ylabel("Error")
legend()


%% LTEs of XB1 and XB2 as a function of step size (h) with fit lines

% compute fit lines
filter_params.min_xval = 0.1;
filter_params.max_xval = 1.5;
[p1, k1] = loglog_fit(h_range, LTE_XB1, filter_params);
[p2, k2] = loglog_fit(h_range, LTE_XB2, filter_params);

figure(2)
% plot data
loglog(h_range, LTE_XB1, ".", MarkerSize=12, Color=C(1), DisplayName=(method_name + ": XB1")); hold on;
loglog(h_range, LTE_XB2, ".", MarkerSize=12, Color=C(2), DisplayName=(method_name + ": XB2"));

% plot fit lines
plot_name = ['XB1 Reg. line p=', num2str(round(p1, 2))];
loglog(h_range, k1*h_range.^p1, "-", Color=C(1), DisplayName=plot_name)
plot_name = ['XB2 Reg. line p=', num2str(round(p2, 2))];
loglog(h_range, k2*h_range.^p2, "-", Color=C(2), DisplayName=plot_name)

% label plot
title("Local Truncation Error vs. Step Size")
xlabel("Step Size"); ylabel("Error")
legend()


%% Plot the LTEs of XB1 and XB2 as a function of their difference, |XB1 − XB2|. 
% How do the errors scale with |XB1 − XB2|? Is |XB1 − XB2| a good proxy for the error 
% (or at least, does it provide a good upper bound on the error?).

figure(3);
loglog(XB_diff, LTE_XB1, ".", MarkerSize=12, DisplayName=(method_name + ": XB1")); hold on;
loglog(XB_diff, LTE_XB2, ".", MarkerSize=12, DisplayName=(method_name + ": XB2"));

% label plot
title("Local Truncation Errors as a function of |XB1 − XB2|")
xlabel("Step Size"); ylabel("Error")
legend()

%% Helper function to get local trunctation error list

function err_list = get_local_tr_err(rate_fcn, XA, V0, t_ref, B_row, h_range, BT_struct, orbit_params)

    % make container
    % err_list = zeros(length(h_range), length(XA));
    err_list = zeros(length(h_range), 1);

    for i = 1:length(h_range)
        
        if B_row == 1
            [X_approx, ~, ~] = explicit_RK_step_embedded(rate_fcn, t_ref, XA, h_range(i), BT_struct);
        else % B row = 2
            [~, X_approx, ~] = explicit_RK_step_embedded(rate_fcn, t_ref, XA, h_range(i), BT_struct);
        end
        
        % compute val at next step
        t = t_ref + h_range(i);
        X_exact = compute_planetary_motion(t, V0, orbit_params);
        
        % store error
        % err_list(i, :) = norm(X_exact - X_approx);
        err_list(i) = norm(X_exact - X_approx);
    end

end


%% Helper function that constructs Butcher Tableaus

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
