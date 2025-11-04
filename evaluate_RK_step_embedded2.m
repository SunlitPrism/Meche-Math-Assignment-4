function [XB1,XB2, difference] = evaluate_RK_step_embedded2(t_ref, h)
    % Testing using DormandPrince embedded RK method for adaptive step
    % sizes
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    rate_func01 = @(t, X) -5*X + 5*cos(t) - sin(t); 
    XA = 1; t_ref = 0.331

    [XB1, XB2, num_evals] = RK_step_embedded2(rate_func01(t_ref, XA),t_ref,XA,h,DormandPrince);
    difference = abs(XB1 - XB2);
    num_evals
end