% This function computes the value of X at the next time step
<<<<<<< Updated upstream
% for any arbitrary embedded RK method
=======
% for any arbitrary Runge KuttA method
>>>>>>> Stashed changes
%
% INPUTS:
%   rate_func_in: the function used to compute dXdt. rate_func_in will
%           have the form: dXdt = rate_func_in(t,X) (t is before X)
%   t: the value of time at the current step
%   XA: the value of X(t)
%   h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%   BT_struct: a struct that contains the Butcher tableau
<<<<<<< Updated upstream
%       BT_struct.A: matrix of a_{ij} values
%       BT_struct.B: vector of b_i values (2-rows for embedded methods)
%       BT_struct.C: vector of c_i values
% OUTPUTS:
%   XB1: the approximate value for X(t+h) using the first row of the Tableau
%   XB2: the approximate value for X(t+h) using the second row of the Tableau
%   num_evals: A count of the number of times that you called
%   rate_func_in when computing the next step

function [XB1, XB2, num_evals] = explicit_RK_step_embedded(rate_func_in, t, XA, h, BT_struct)
=======
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
% OUTPUTS:
%   XB: the approximate value for X(t+h) (the next step)
%           formula depends on the integration method used
%   num_evals: A count of the number of times that you called
%   rate_func_in when computing the next step

% takes into account an embedded step (2 rows in BT_struct.B)
function [XB, num_evals] = explicit_RK_step_embedded(rate_func_in, t, XA, h, BT_struct)
>>>>>>> Stashed changes
    
    num_evals = 0;
    K = zeros(length(XA), length(BT_struct.B));

    for n = 1:length(BT_struct.B)
        t_temp = t + BT_struct.C(n)*h;
        X_temp = XA + h*(K*BT_struct.A(n,:)')';
        K(:, n) = rate_func_in(t_temp, X_temp);

        num_evals = num_evals+1;
    end

    XB1 = XA + h*(K*BT_struct.B(1,:)')';
    XB2 = XA + h*(K*BT_struct.B(2,:)')';
<<<<<<< Updated upstream
=======
    XB = [XB1; XB2]
>>>>>>> Stashed changes

end