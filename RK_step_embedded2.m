%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded2(rate_func_in,t,XA,h,BT_struct)
    % B1 = BT_struct.B(:,1);
    % B2 = BT_struct.B(:,2);
    % k_B1 = zeros(length(XA), length(B1));
    % k_B2 = zeros(length(XA), length(B2));
    % 
    % num_evals = 0;
    % 
    % sum_B1 = 0;
    % sum_B2 = 0;
    % for i = 1:length(B1)
    %     t_input = t + BT_struct.C(i) * h;
    %     X_input = 0;       
    %     for j = 1:i-1
    %         X_input = X_input + BT_struct.A(i, j) * k_B1(:,j );
    %     end
    %     X_input = XA + h*X_input;
    % 
    %     k_B1(:, i) = rate_func_in(t_input, X_input);
    % 
    %     num_evals = num_evals + 1;
    % 
    %     sum_B1 = sum_B1 + B1(i) * k_B1(:, i);
    % end
    % 
    % for i = 1:length(B2)
    %     t_input = t + BT_struct.C(i) * h;
    %     X_input = 0;       
    %     for j = 1:i-1
    %         X_input = X_input + BT_struct.A(i, j) * k_B2(:,j );
    %     end
    %     X_input = XA + h*X_input;
    % 
    %     k_B2(:, i) = rate_func_in(t_input, X_input);
    % 
    %     num_evals = num_evals + 1;
    % 
    %     sum_B2 = sum_B2 + B2(i) * k_B2(:, i);
    % end
    % 
    % XB1 = XA + h*sum_B1;
    % XB2 = XA + h*sum_B2;
    k = zeros(length(XA), length(BT_struct.B));
    sum = 0;
    num_evals = 0;

    for i = 1:length(BT_struct.B)
        
        t_input = t + BT_struct.C(i) * h;
        X_input = 0;       
        for j = 1:i-1
            X_input = X_input + BT_struct.A(i, j) * k(:,j );
        end
        X_input = XA + h*X_input;

        k(:, i) = rate_func_in(t_input, X_input);

        num_evals = num_evals + 1;

        % sum = sum + BT_struct.B(1,i) * k(:, i);
    end
 
    XB1 = XA + h*k*BT_struct.B(1,:)';
    XB2 = XA + h*k*BT_struct.B(2,:)';

end