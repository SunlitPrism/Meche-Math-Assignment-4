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
function [t_list, X_list, h_avg, h_next_list, num_evals, step_fail_rate] = explicit_RK_variable_step_integration ...
    (rate_func_in, tspan, X0, h_ref, BT_struct, p, error_desired)
    
    N = ceil(diff(tspan)/h_ref)+1; % add 1 bc of inclusivity in linspace range

    t_list = linspace(tspan(1), tspan(2), N);
    X_list = zeros(N, length(X0));

    h_init = diff(t_list(1:2));  
    h_next_list = h_init;
    
    total_evals = 0; % counter for num times rate func is called

    % start with inital condition, calc ans store future vals
    X_list(1, :) = X0;
    redo = false;
    failed_step = 0;
    all_step = 0;

    for i = 1:N-1

        % select t value
        t = t_list(i);

        if redo == true
            
            disp("redo iter")
            failed_step = failed_step+1;
            all_step = all_step+1;
            
            % find x_(i+1), update x_i (XA here)
            [XB, num_evals, h_next, redo] = explicit_RK_variable_step ... 
                        (rate_func_in, t, X_list(i,:), h_next_list(i), BT_struct, p, error_desired);
            
            % update and store h value
            h_next_list(i+1) = h_next;
            total_evals = total_evals + num_evals;

        elseif redo == false
            
            % find x_(i+1), update x_i (XA here)
            all_step = all_step+1;
            [XB, num_evals, h_next, redo] = explicit_RK_variable_step ...
                        (rate_func_in, t, X_list(i,:), h_next_list(i), BT_struct, p, error_desired);
            
            % update and store h value
            h_next_list(i+1) = h_next;

            % update and store vals
            X_list(i+1, :) = XB(1,:); % storing just XB1 
            total_evals = total_evals + num_evals;

        else
            disp('error')
        end
    end

    % compute step failure rate and the average step size
    step_fail_rate = failed_step/all_step;
    h_avg = mean(h_next_list);

end