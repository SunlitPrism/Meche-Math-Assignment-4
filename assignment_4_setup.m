function assignment_4_setup()
    t = 0.5;
    x_p = 0;
    y_p = 1;
    % dxdt_p =
    % dydt_p =
    % V = [x_p; y_p; dxdt_p; dydt_p];
    % orbit_params = struct();
    % orbit_params.m_sun = 
    % orbit_params.m_planet =
    % orbit_params.G = 6.6743e-11; %m3kg-1 s-2

    rate_func_in = @rate_func02; %gravity_rate_func;

    XA = [x_p,y_p];
    h = 0.1;


    BT_struct = struct();
    BT_struct.A = [0,0; 0.5,0]; %Matrix 
    BT_struct.B = [0;0.5]; %Vector
    BT_struct.C = [0;1]; %Vector

    disp(BT_struct)
    
    [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct);

    disp(XB)
    num_evals(XB)
end