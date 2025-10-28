%%

% choose three or more different explicit Runge-Kutta methods to test

% Ideally, there will be a few different orders of truncation error that 
% are represented from among the methods that you chose.

% Perform a set of experiments comparing the behavior of the different method

%  plot comparing the approximate solution for each method to the true solution for a given time step-length.
% • A log log plot showing how the local truncation error (and the difference |X(t + h) − X(t)|) scales with the
% step-size, h (with fit-lines).
% • A log log plot showing how the global truncation error scales with the step-size, h (with fit-lines).
% • A log log plot showing how the global truncation error scales with num evals (with fit-lines).

%% Initialize several Butcher Tableaus 

method_list = ["Dormand Prince", "Fehlberg", "Heun Euler", "Fehlberg RK1", "Bogacki"];
%selection = logical([1, 1, 1, 1, 1]);
selection = logical([0, 0, 0, 1, 1]);
selected_methods = method_list(selection); 


num_methods = sum(selection);
BT_list = cell(1, num_methods);

% populate BT_list for the selected methods
for i = 1:num_methods
    BT_list{i} = get_BT(selected_methods(i));
end


%% Initialize system

% orbit params
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

% initial state
x0 = 8; y0 = 1; 
<<<<<<< Updated upstream
dxdt0 = 0.1; dydt0=0;
dxdt0 = 0; dydt0=0;
=======
>>>>>>> Stashed changes
dxdt0 = 0; dydt0=0;
V0 = [x0, y0, dxdt0, dydt0];

%% Compare approx. soln to true soln for various methods
% construct a subplot for each method 

<<<<<<< Updated upstream
clf; h_ref = 0.01; tspan = [0, 30]; 
clf; h_ref = 0.2; tspan = [0, 30]; 
=======
>>>>>>> Stashed changes
clf; h_ref = 0.2; tspan = [0, 30]; 
t_range = linspace(tspan(1), tspan(2), ceil(diff(tspan)/h_ref)+1);
h = diff([t_range(1), t_range(2)]);

% compute true soln
V_list = compute_planetary_motion(t_range, V0, orbit_params);

% anonymize func to pass as input later
my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);

<<<<<<< Updated upstream
figure(1);
for j = 1:num_methods

    subplot(num_methods, 1, j)
    
    % plot "true" solution
    plot(t_range, V_list(:,1), DisplayName="True soln x"); hold on;
    % plot(t_range, V_list(:,2), DisplayName="True soln y");

end

=======
>>>>>>> Stashed changes
subplot(2, 1, 1)
    
% plot "true" solution
plot(t_range, V_list, DisplayName="True soln"); 

%figure(1);
for j = 1:num_methods

    subplot(num_methods, 1, j)
    method_name1 = selected_methods(j)
    
    % plot "true" solution
    plot(t_range, V_list, DisplayName="True soln"); hold on;

    % solve with numerical method
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration ...
                            (my_rate, tspan, V0, h_ref, BT_list{j});

    % plot numerical
    method_name = selected_methods(j);
<<<<<<< Updated upstream
    plot(t_list, X_list(:,1), ".-", DisplayName=(method_name + " x")); 
    % plot(t_list, X_list(:,2), ".-", DisplayName=(method_name + " y")); 

    % label graph
    method_name = selected_methods(j);
    title("True Soln vs. " + method_name + " Approximation")
    xlabel("Time"); ylabel("X(t)")
    legend()

end

sgtitle("Comparison of Runge-Kutta methods (\Deltat=" + num2str(h) + ") to True Soln")

=======
    plot(t_range, V_list, ".-", DisplayName=method_name); 

    % label graph
    title("True Soln vs. ", method_name, " Approx.")
    xlabel("Time"); ylabel("X(t)")

end

% label plot
sgtitle("Comparison of Runge-Kutta methods (\Delta t=", h, ") to True Soln")
>>>>>>> Stashed changes

%% Local truncation error

% set h range
% logspace

% compute |X(t + h) − X(t)|

% for each method
% compute error
% plot

% label plot

disp('temp stop pt')


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


