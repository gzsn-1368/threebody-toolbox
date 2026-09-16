function [Xdot, jac] = crtbp(~, X, u)
%CRTBP Compute the state derivative and Jacobian of the circular restricted
%   three-body problem (CRTBP) in the dimensionless Earth--Moon rotating
%   coordinate system.
%
%   XDOT = CRTBP(T, X) computes the state derivative at time T for the
%   normalized Earth--Moon mass ratio MU. The time argument is accepted for
%   ODE-solver compatibility and is not used by the circular model.
%
%   [XDOT, JAC] = CRTBP(T, X) also returns the 6-by-6 state Jacobian JAC,
%   which can be used to linearize the dynamics or integrate a state
%   transition matrix.
%
%   [XDOT, JAC] = CRTBP(T, X, U) specifies an optional three-dimensional
%   control acceleration U. If U is omitted, it defaults to zero.
%
%   If X has 42 elements, the first six elements are the state
%   [x; y; z; dx; dy; dz] and the remaining 36 elements are a column-wise
%   6-by-6 state transition matrix PHI. The variational equation
%   dPHI/dt = JAC*PHI is then appended to XDOT.
%
%   Inputs:
%     T   - Scalar time argument, accepted for ODE-solver compatibility.
%     X   - State vector with 6 or 42 elements.
%     U   - Optional 3-vector of dimensionless control acceleration.
%
%   Outputs:
%     XDOT - Column vector of state and, for a 42-element input, variational
%            equation derivatives.
%     JAC  - 6-by-6 state Jacobian.
%
%
%   Copyright 2025 Beijing Jiaotong University.

    if nargin < 3
        u = zeros(3, 1);
    end

    persistent mu;
    if isempty(mu)
        mu = constants('MassUnitEM');
    end
    
    nX = length(X);
    if nX ~= 6 && nX ~= 42
        error('crtbp:InvalidInput', 'X must contain 6 or 42 elements.');
    end
    if length(u) ~= 3
        error('crtbp:InvalidInput', 'u must contain 3 elements.');
    end

    % Equations of motion.
    x = X(1);
    y = X(2);
    z = X(3);
    dx = X(4);
    dy = X(5);
    dz = X(6);

    % Control acceleration.
    ux = u(1);
    uy = u(2);
    uz = u(3);

    r1cubic = ((mu + x).^2 + y.^2 + z.^2).^(3/2);
    r2cubic = ((mu + x - 1).^2 + y.^2 + z.^2).^(3/2);

    Xdot = [
        dx;
        dy;
        dz;
        2*dy + x - (x+mu)*(1-mu)/r1cubic - (x+mu-1)*mu/r2cubic + ux;
        -2*dx + y - y*(1-mu)/r1cubic - y*mu/r2cubic + uy;
        -z*(1-mu)/r1cubic - z*mu/r2cubic + uz
        ];

    % First-order variational equations.
    if nX == 42 || nargout > 1
        r1power5 = ((mu + x)^2 + y^2 + z^2)^(5/2);
        r2power5 = ((mu + x - 1)^2 + y^2 + z^2)^(5/2);
        three_mu = 3 * mu;

        dxdx = 1 + 3*(1-mu)*(x+mu)^2/r1power5 + three_mu*(1-mu-x)^2/r2power5 - ...
            (1-mu)/r1cubic - mu/r2cubic;
        dydy = 1 + 3*(1-mu)*y^2/r1power5 + three_mu*y^2/r2power5 - ...
            (1-mu)/r1cubic - mu/r2cubic;
        dzdz = 3*(1-mu)*z^2/r1power5 + three_mu*z^2/r2power5 - (1-mu)/r1cubic - ...
            mu/r2cubic;
        dxdy = 3*(1-mu)*(x+mu)*y/r1power5 - three_mu * (1-mu-x)*y/r2power5;
        dydx = dxdy;
        dxdz = 3*(1-mu)*(x+mu)*z/r1power5 - three_mu * (1-mu-x)*z/r2power5;
        dzdx = dxdz;
        dydz = 3*(1-mu)*y*z/r1power5 + three_mu*y*z/r2power5;
        dzdy = dydz;

        jacobiLowerLeft = [
            dxdx, dxdy, dxdz;
            dydx, dydy, dydz;
            dzdx, dzdy, dzdz
            ];
        jacobian = [
            zeros(3,3), eye(3); 
            jacobiLowerLeft, [0, 2, 0; -2, 0, 0; 0, 0, 0]
            ];

        if nX == 42
            stateTransitionMatrix = jacobian * reshape(X(7:42), 6, 6);
            Xdot = [Xdot; reshape(stateTransitionMatrix, 36, 1)];
        end
        if nargout > 1
            jac = jacobian;
        end
    end
end
