function assignment_4_setup_avery_saving_backup()

    %Establish Constants/variables
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

    %Compute Motion
    t_range = linspace(t0,tf,t_iter);
    V_list = compute_planetary_motion(t_range, V0, orbit_params);

    %Set up Rate Function
    rate_func_in = @(t_in,V_in) gravity_rate_func(t_in,V_in,orbit_params);

    h_ref = 0.01;

    BT_struct = struct();
    BT_struct.A = [0,0; 0.5,0]; %Matrix 
    BT_struct.B = [0;0.5]; %Vector
    BT_struct.C = [0;1]; %Vector

    %Explicit RK function
    [t_list, X_list, h_avg, num_evals]= explicit_RK_fixed_step_integration(rate_func_in,tspan,V0,h_ref,BT_struct);

    %Plot Velocity List vs Time and Postion List vs Time
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

    %Errors vs H size Calculations Global
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

    %Fit and Regression Global
    filter_params = struct();
    filter_params.min_yval = 1e-10;
    filter_params.max_yval = 1;
    
    [p1,k1] = loglog_fit(h_avg_list,tr_error_list,filter_params);
    [p2,k2] = loglog_fit(h_avg_list,tr_error_list,filter_params);

    p1 = abs(p1);
    p2 = abs(p2);

    %Plot Global Error as Function of Time
    figure(2);
    loglog(h_avg_list,tr_error_list,'ro','MarkerFaceColor','r') %Global as func of time
    loglog(h_avg_list,k1*h_avg_list.^p1,'k') %Global as func of time

    %Plot Global Error as Function of Evaluations
    figure(3);
    loglog(num_evals_list,tr_error_list,'bo','MarkerFaceColor','b') %Global as #evals
    loglog(num_evals_list,k2*num_evals_list.^p2,'k') %Global as #evals

     %Errors vs H size Calculations Local
    n_samples = 60;
    h_ref_list = logspace(-3,-1,n_samples);

    abs_diff_list = [];
    tr_error_list1 = [];
    tr_error_list2 = [];

    for n = 1:length(h_ref_list)
       h_ref = h_ref_list(n);
       V_list = compute_planetary_motion(tspan(1)+h_ref,V0,orbit_params);

       [XB1,XB2,~]= explicit_RK_step_embedded(rate_func_in,tspan,V0,h_ref,BT_struct);
       abs_diff_list(n) = norm(V_list - V0);
        tr_error_list1(n) = norm(XB1-V_list);
        tr_error_list2(n) = norm(XB2-V_list);
    end
     %Fit and Regression Local
    filter_params = struct();
    filter_params.min_yval = 1e-10;
    filter_params.max_yval = 1;
    
    [p1,k1] = loglog_fit(h_avg_list,tr_error_list1,filter_params);
    [p2,k2] = loglog_fit(h_avg_list,tr_error_list2,filter_params);

    p1 = abs(p1);
    p2 = abs(p2);

    %Plot Global Error as Function of Time
    figure(2);
    loglog(h_avg_list,tr_error_list,'ro','MarkerFaceColor','r') %Global as func of time
    loglog(h_avg_list,k1*h_avg_list.^p1,'k') %Global as func of time

    %Plot Global Error as Function of Evaluations
    figure(3);
    loglog(num_evals_list,tr_error_list,'bo','MarkerFaceColor','b') %Global as #evals
    loglog(num_evals_list,k2*num_evals_list.^p2,'k') %Global as #evals
end