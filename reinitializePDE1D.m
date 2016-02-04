function phi = reinitializePDE1D(x,psi,delay)
%   reinitializes the level set psi(x) into a signed-distance function
%   by solving
%
%       phi_tau(x,tau) - sgn(psi(x))*(1 - |phi_x(x,tau)|) = 0
%       phi(x,0) = psi(x)
%
%   until numerical steady-state
%
%   inputs:
%       x - vector containing x_i
%       psi - vector containing psi(x_i,t)
%       delay - animation delay (optional)
%
%   outputs:
%       phi - vector containing phi(x_i,t)

    if (nargin == 2)
        delay = 0;
    end
    
    dx = x(2) - x(1);
    
    % set algorithm parameters
    % (adopted from Sussman, Smereka, Osher, JCP 1994)
    stopping_tolerance = dx^2; % dtau multiplied in while condition
    smoothing_epsilon = dx;
    CFL = 0.9;
    
    % advance phi in tau-time until stopping criterion is met
    % equation first rewritten as phi_tau + (sgn*|phi_x|/phi_x)*phi_x = sgn
    sgn_vals = sgn(psi,smoothing_epsilon);
    phi_old = psi + 2*stopping_tolerance;
    phi_new = psi;
    dtau = 0;
    s = size(x);
    s(s > 1) = s(s > 1) + 1;
    u = zeros(s);
    while (norm(phi_new-phi_old,1)/length(x) > dtau*stopping_tolerance)
        phi_old = phi_new;
        
        % compute advection velocity at cell edges using extrapolation BCs
        phi_x = (phi_old(2:end) - phi_old(1:end-1))/dx;    
        u(2:end-1) = 0.5*(sgn_vals(1:end-1)+sgn_vals(2:end)) ...
                        .*abs(phi_x)./phi_x;
        u(1) = 2.0*u(2) - u(3);
        u(end) = 2.0*u(end-1) - u(end-2);
        dtau = CFL*dx/max(abs(u)); % adaptive timestep
        
        phi_new = advect1D(x,phi_old,u,dtau,sgn_vals);
        if (delay > 0)
            plot(x,phi_new,'-b',x,psi,'--r')
            pause(delay);
        end
    end
       
    phi = phi_new;
end

function val = sgn(z,epsilon)
%   function to evaluate the smoothed signum function
%
%   sgn(z) = z / sqrt(z^2 + epsilon)

    val = z ./ sqrt(z.^2 + epsilon);
    
end