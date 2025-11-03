%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will have the form: dXdt = rate_func_in(t,X) (t is before X)
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
function [t_list, X_list, h_avg, h_list, total_evals, step_fail_rate] = explicit_RK_variable_step_integration ...
    (rate_func_in, tspan, X0, h_ref, BT_struct, p, error_desired)

    t_list = [];
    X_list = [];
    h_list = [];

    % initialize lists
    t_list(1) = tspan(1);
    X_list(:, 1) = X0;

    h = h_ref; % set first value
    h_list(1) = h; % store first value

    % initialize counters
    total_evals = 0; 
    failed_step = 0;
    all_step = 0;
    i = 1;

    while t_list(end) < tspan(2)

        % find x_(t + h)
        [XB, num_evals, h_next, redo] = explicit_RK_variable_step ... 
                        (rate_func_in, t_list(i), X_list(:, i), h, BT_struct, p, error_desired);

        % update counters
        all_step = all_step+1;
        total_evals = total_evals + num_evals;
        
        if redo == false

            % update and store vals
            X_list(:, i+1) = XB; % storing just XB1 
            t_list(i+1) = t_list(i) + h;
            
            h = h_next;
            h_list(i+1) = h;
            
            % increment step
            i = i+1;
        else

            h = h_next; % update h
            failed_step = failed_step + 1; % update counter

        end
    end

    % compute step failure rate and the average step size
    step_fail_rate = failed_step/all_step;
    h_avg = mean(h_list);

end