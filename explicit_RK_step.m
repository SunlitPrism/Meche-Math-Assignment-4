% This function computes the value of X at the next time step
% for any arbitrary Runge KuttA method
%
% INPUTS:
%   rate_func_in: the function used to compute dXdt. rate_func_in will
%           have the form: dXdt = rate_func_in(t,X) (t is before X)
%   t: the value of time at the current step
%   XA: the value of X(t)
%   h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%   BT_struct: a struct that contains the Butcher tableau
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
% OUTPUTS:
%   XB: the approximate value for X(t+h) (the next step)
%           formula depends on the integration method used
%   num_evals: A count of the number of times that you called
%   rate_func_in when computing the next step

function [XB, num_evals] = explicit_RK_step(rate_func_in, t, XA, h, BT_struct)
    
    num_evals = 0;
    K = zeros(length(XA), length(BT_struct.B));
    % i made ours wrong, its like 5x1 instead of 1x5 for example

    for n = 1:length(BT_struct.B)
        t_temp = t + BT_struct.C(n)*h;
        X_temp = XA + h*(K*BT_struct.A(n,:)')';
        K(:, n) = rate_func_in(t_temp, X_temp);

        num_evals = num_evals+1;
    end

    XB = XA + h*(K*BT_struct.B')';

    % % initialize values and containers
    % s = size(BT_struct, 1);
    % k_vals = zeros(s, length(XA));
    % k_vals(1, :) = rate_func_in(t, XA); % dxdt at initial point
    % 
    % % compute all k values
    % for i = 2:s
    % 
    %     % index the constants
    %     ci = BT_struct.C(i);
    %     a_vals = BT_struct.A(i, 1:i-1);
    %     k_in = k_vals(1:i-1, :);
    % 
    %     % calc and store next k value
    %     k_vals(i, :) = rate_func_in(t+ci, XA + a_vals.*k_in);
    % 
    % end
    % 
    % % calculate E(i=1, s) for b_i*k_i
    % biki = BT_struct.B .* k_vals;
    % 
    % % compute next step
    % XB = XA + h*biki;
    % num_evals = s;

end