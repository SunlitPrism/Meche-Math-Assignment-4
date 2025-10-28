%% in class following along with orion

function class_10_21()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    x0 = 8; y0 = 0; 
    dxdt0 = 0;
    dydt0=0;

    V0 = [x0, y0, dxdt0, dydt0];
    t_range = linspace(0,  30, 100);
    V_list = compute_planetary_motion(t_range, V0, orbit_params);

    my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);


    BT_struct = []; %% grab info from assignment page
    h_ref = 0.1;
    [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration ...
                            (my_rate, tspan, V0, h_ref, BT_struct);

    subplot(2, 1, 1);
    hold on;
    plot(t_range, V_list(:,1), "k", LineWidth=2)
    plot(t_range, V_list(:,2), "b", LineWidth=2)

    plot(t_list, X_list(:,1), "r--", LineWidth=2)
    plot(t_list, X_list(:,2), "r--", LineWidth=2)

    subplot(2, 1, 2);
    hold on;
    plot(t_range, V_list(:,3), "k", LineWidth=2)
    plot(t_range, V_list(:,4), "b", LineWidth=2)

    plot(t_list, X_list(:,3), "r--", LineWidth=2)
    plot(t_list, X_list(:,4), "r--", LineWidth=2)

    num_samples = 30;
    h_ref_list = logspace(-5, -1, num_samples);

    num_evals_list = [];
    h_avg_list = [];
    tr_error_list = [];

    for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);

        [~, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration ...
                            (my_rate, tspan, V0, h_ref, BT_struct);

        tr_error = norm(X_list(end, :)-V_list(end,:));
        tr_error_list(n) = tr_error;
        h_avg_list(n) = h_avg;
        num_evals_list(n) = num_evals;

    end

    filter_params = struct(); filter_params.min_

    [p1, k1] = loglog_fit(h_avg_list, tr_error_list, filter_params);
    [p2, k2] = loglog_fit(num_evals_list, tr_error_list, filter_params);

    % display
    abs(p1) 
    abs(p2)

    % plot global tr vs step size AND num evals
    figure(2);
    loglog(h_avg_list, tr_error_list, "ro", MarkerFaceColor="r")
    hold on;
    loglog(h_avg_list, k1*h_avg_list^p1, "k")

    figure(3);
    loglog(num_evals_list, tr_error_list, "ro", MarkerFaceColor="r")
    loglog(num_evals_list, k2*num_evals_list^p2, "k")

    % embedded experiment
    % green line
    figure(2); 
    loglog(h_ref_list, abs_diff_list, "ro", MarkerFaceColor="r")
    hold on;
    loglog(h_ref_list, tr_error_list1, "go", MarkerFaceColor="g")
    loglog(h_ref_list, tr_error_list2, "bo", MarkerFaceColor="b")

    % compute fits again
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    filter_params = struct(); % define again, ymX_list(1,:) = X0';in = 10^-13, ymax_10^-6
=======
    filter_params = struct(); % define again, ymin = 10^-13, ymax_10^-6
>>>>>>> Stashed changes
=======
    filter_params = struct(); % define again, ymin = 10^-13, ymax_10^-6
>>>>>>> Stashed changes
    [p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
    [p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);

    % plot reg lines
    loglog(h_ref_list, k1*h_ref_list^p1, "k")
    loglog(h_ref_list, k2*h_ref_list^p2, "k")


end