%% graphs for write up (testing section)

% choose three or more different explicit Runge-Kutta methods to test

% Ideally, there will be a few different orders of truncation error that 
% are represented from among the methods that you chose.

% Perform a set of experiments comparing the behavior of the different method

%  plot comparing the approximate solution for each method to the true solution for a given time step-length.
% • A log log plot showing how the local truncation error (and the difference |X(t + h) − X(t)|) scales with the
% step-size, h (with fit-lines).
% • A log log plot showing how the global truncation error scales with the step-size, h (with fit-lines).
% • A log log plot showing how the global truncation error scales with num evals (with fit-lines).

%% Housekeeping

clc; clear;

%% Initialize several Butcher Tableaus 

method_list = ["Dormand Prince", "Fehlberg", "Heun Euler", "Fehlberg RK1", "Bogacki"];
selection = logical([1, 0, 1, 0, 1]);
p_method = [5,4; 5,4; 2,1; 2,1; 3,2;]';
selected_methods = method_list(selection); 

num_methods = sum(selection);
BT_list = cell(1, num_methods);

% populate BT_list for the selected methods
for i = 1:num_methods
    BT_list{i} = get_BT(selected_methods(i));
end


%% Initialize system

% first order system driven by a sinusoid
rate_func01 = @(t, X) -5*X + 5*cos(t) - sin(t); V0 = 1;
solution01 = @(t) cos(t);

my_rate = rate_func01;

%% Compare approx. soln to true soln for various methods
% construct a plot for each method 

C = rgb2hex(orderedcolors("gem"));

% h_ref = 0.3; tspan = [0, 8]; 
% t_range = linspace(tspan(1), tspan(2), ceil(diff(tspan)/0.1)+1);
% 
% % compute true soln
% V_list = solution01(t_range);
% my_rate = rate_func01;
% % my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);
% 
% for j = 1:num_methods
% 
%     % initialize figure
%     figure(j);
% 
%     % solve with numerical method
%     [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration ...
%                             (my_rate, tspan, V0, h_ref, BT_list{j});
% 
%     % plot "true" solution
%     plot(t_range, V_list, "-", Color=C(3), LineWidth=2, DisplayName="Analytical soln");
%     hold on;
% 
%     % plot numerical
%     method_name = selected_methods(j);
%     plot(t_list, X_list, ".-", DisplayName=(method_name), Color=C(1), MarkerSize=12); 
%     hold off;
% 
%     % label graph
%     h_str = num2str(round(h_avg, 3));
%     title(method_name + " Approx. (\Deltat=" + h_str + ") vs. Analytical Solution")
%     xlabel("Time"); ylabel("X(t)"); ylim([-1.2, 1.2]);
%     legend()
% 
% end


%% Local truncation error

% set h range
h_len = 100; t_ref = 0.331;
h_range = logspace(-6, 2, h_len);

% called above, here for reference
% my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);
X0 = solution01(t_ref);

% compute |X(t + h) − X(t)|
ref_X_list = zeros(h_len, length(X0));

% get ref line data
for i = 1:h_len
    t_vals = [t_ref, t_ref + h_range(i)];
    V_list = solution01(t_vals);
    ref_X_list(i,:) = abs(diff(V_list));
end

filter_params.max_xval = 0.9;

clf
figure(1);
% plot ref line 
loglog(h_range, ref_X_list(:,1), ".", Color=C(1), DisplayName="Ref. line"); hold on;

% calc & plot loglog regression fit line
[p_ref, k_ref] = loglog_fit(h_range, ref_X_list, filter_params);
plot_name = ['Ref. line p=', num2str(round(p_ref, 2))];
loglog(h_range, k_ref*h_range.^p_ref, "-", Color=C(1), DisplayName=plot_name)

% for each RK method compute local tr error

for j = 1:num_methods

    % get local truncation error for jth method
    Xlist = get_local_tr_err(my_rate, solution01, X0, t_ref, h_range, BT_list{j}, 1);

    % plot results
    method_name = selected_methods(j);
    loglog(h_range, Xlist, ".", Color=C(j+1), DisplayName=method_name)

    % plot loglog regression line
    [p, k] = loglog_fit(h_range, Xlist, filter_params);
    plot_name = ['Reg. line p=', num2str(round(p, 2))];
    loglog(h_range, k*h_range.^p, "-", Color=C(j+1), DisplayName=plot_name)

end

% label plot
title("Local Truncation Error for RK Methods")
xlabel("Step size"); ylabel("Error"); 
legend()



%% Helper function to get local trunctation error list

function err_list = get_local_tr_err(rate_fcn, solver, XA, t_ref, h_range, BT_struct, B_row)

    % make container
    err_list = zeros(length(h_range), length(XA));

    for i = 1:length(h_range)

        % compute val at next step
        t = t_ref + h_range(i);

        if B_row == 1
            [X_approx, X1, ~] = RK_step_embedded2(rate_fcn, t, XA, h_range(i), BT_struct);
        else
            [X1, X_approx, ~] = RK_step_embedded2(rate_fcn, t, XA, h_range(i), BT_struct);
        end

        % debugging
        if X_approx == X1
            disp("XB1 = XB2")
        end

        X_exact = solver(t);
        
        % store error
        err_list(i, :) = abs(X_exact - X_approx);
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


