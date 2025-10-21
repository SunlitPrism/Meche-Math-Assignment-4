function assignment_4_setup_avery_saving_backup()
    t0 = 0.5;
    tf = 6;
    tspan = [t0,tf];
    t_iter = 100;
    x_p = 8;
    y_p = 0;
    dxdt_p = 0;
    dydt_p = 1.5;
    V0 = [x_p; y_p; dxdt_p; dydt_p];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    t_range = linspace(t0,tf,t_iter);
    V_list = compute_planetary_motion(t_range, V0, orbit_params);

    rate_func_in = @(t_in,V_in) gravity_rate_func(t_in,V_in,orbit_params);

    h_ref = 0.01;

    BT_struct = struct();
    BT_struct.A = [0,0; 0.5,0]; %Matrix 
    BT_struct.B = [0;0.5]; %Vector
    BT_struct.C = [0;1]; %Vector
    
    [t_list, X_list, h_avg, num_evals]= explicit_RK_fixed_step_integration(rate_func_in,tspan,V0,h_ref,BT_struct);

    figure(1);
    subplot(2,1,1);
    hold on;
    plot(t_range,V_list(:,1));
    plot(t_range,V_list(:,2));

    plot(t_list,X_list(:,1));
    plot(t_list,X_list(:,2));
    hold off;

    subplot(2,1,2)
    hold on;
    plot(t_range,V_list(:,3));
    plot(t_range,V_list(:,4));

    plot(t_list,X_list(:,3));
    plot(t_list,X_list(:,4));
    hold off;

    h_ref_list = logspace(-3,-1,30);

    num_evals_list = [];
    h_avg_list = [];
    tr_error_list = [];

    for n = 1:length(h_ref_list)
       h_ref = h_ref_list(n);

       [t_list, X_list, h_avg, num_evals]= explicit_RK_fixed_step_integration(rate_func_in,tspan,V0,h_ref,BT_struct);

       tr_error = norm(X_list(end,:)-V_list(:,end));
       tr_error_list(n) = tr_error;
       h_avg_list(n) = h_avg;
       num_evals_list(n) = num_evals;
    end

    filter_params = struct();
    % filter_params.min_xval;
    % filter_params.max_xval;
    filter_params.min_yval = 1e-10;
    filter_params.max_yval = 1;

    [p1,k1] = loglog_fit(h_avg_list,tr_error_list,filter_params);
    [p2,k2] = loglog_fit(h_avg_list,tr_error_list,filter_params);

    p1 = abs(p1);
    p2 = abs(p2);

    figure(2);
    loglog(h_avg_list,tr_error_list,'ro','MarkerFaceColor','r') %Global as func of time
    loglog(h_avg_list,k1*h_avg_list.^p1,'k') %Global as func of time

    figure(3);
    loglog(num_evals_list,tr_error_list,'bo','MarkerFaceColor','b') %Global as #evals
    loglog(num_evals_list,k2*num_evals_list.^p2,'k') %Global as #evals
end