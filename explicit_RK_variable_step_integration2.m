

%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg,step_failure_rate, num_evals] = explicit_RK_variable_step_integration2(rate_func_in,tspan,X0,h_ref,DormandPrince,p,error_desired)
    t_list = [];
    X_list = [];
    h_vals = [];
    attempt_steps = 0;
    failed_steps = 0;
    num_evals = 0;
    
    t_list(1) = tspan(1);
    X_list(:,1) = X0;
    i = 2;
    h = h_ref;
    h_vals(1) = h_ref;


    while t_list(end) < tspan(end)
        % now using our variable step
        [XB, num_eval_temp, h_next, redo] = explicit_RK_variable_step(rate_func_in, t_list(i-1), X_list(:, i-1), h, DormandPrince, p, error_desired);
        num_evals = num_evals+num_eval_temp;
        
        if redo == false
            X_list(:, i) = XB;
            t_list(i) = t_list(i-1) + h;
            h = min(h_next, tspan(end)-t_list(i)+10^(-15));
            attempt_steps = attempt_steps + 1;
            
            h_vals(i) = h;
            i = i+1;
        else
            h = h_next;
            failed_steps = failed_steps + 1;
            attempt_steps = attempt_steps + 1;
        end

    end

    h_avg = mean(h_vals);
    step_failure_rate = failed_steps/attempt_steps;
end