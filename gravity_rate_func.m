% Rate function describing Newton's law of gravitation
%
% INPUTS:
%   t: the current time
%   V: the vector of the position and velocity of the planet
%   V = [x_p; y_p; dxdt_p; dydt_p]
%   orbit_params: a struct describing the system parameters
%       orbit_params.m_sun is the mass of the sun
%       orbit_params.m_planet is the mass of the planet
%       orbit_params.G is the gravitational constant
% OUTPUTS:
%   dVdt: a column vector describing the time derivative of V:
%   dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]

function dVdt = gravity_rate_func(t, V, orbit_params)
    
    % unpack variables
    G = orbit_params.G;
    ms = orbit_params.m_sun;
    mp = orbit_params.m_planet;
    r = [V(:,1), V(:,2)]; % [x_p, y_p], position of planet

    % Newton's Law of Universal Gravitation
    % m_p * (d^2r/dt^2) = F = - (m_p*m_s*G)/|r^3| * r

    % calc second deriv
    dr2dt2 = -(ms*G/norm(r^3))*r;

    % isolate x and y components
    d2xdt2_p = dr2dt2(1);
    d2ydt2_p = dr2dt2(2);

    % construct dVdt
    dVdt = [V(:,3); V(:,4); d2xdt2_p; d2ydt2_p];

end
