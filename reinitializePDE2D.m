function phi = reinitializePDE2D(x,y,psi,delay)
%   reinitializes the level set psi(x,y) into a signed-distance function
%   by solving
%
%       phi_tau(x,y,tau) - sgn(psi(x,y))*(1 - |grad_phi(x,y,tau)|) = 0
%       phi(x,y,0) = psi(x,y)
%
%   until numerical steady-state
%
%   inputs:
%       x - vector containing x_j
%       y - vector containing y_i
%       psi - vector containing psi(x_j,y_i,t)
%       delay - animation delay (optional)
%
%   outputs:
%       phi - vector containing phi(x_j,y_i,t)

    if (nargin == 2)
        delay = 0;
    elseif (delay > 0)
        close all;
    end
    
    dx = x(1,2) - x(1,1);
    
    % set algorithm parameters
    % (first two are adopted from Sussman, Smereka, Osher, JCP 1994)
    stopping_tolerance = dx^2; % dtau multiplied in while condition
    smoothing_epsilon = dx;
    CFL = 0.9;
    
    % advance phi in tau-time until stopping criterion is met
    % equation first rewritten as phi_tau + 
    % (sgn*grad_phi/|grad_phi|)*dot*grad_phi = sgn
    sgn_vals = sgn(psi,smoothing_epsilon);
    phi_old = psi + 2*stopping_tolerance;
    phi_new = psi;
    dtau = 0;
    phi_x = zeros(size(x));
    phi_y = zeros(size(x));
    while (norm(phi_new-phi_old,1)/numel(x) > dtau*stopping_tolerance)
        phi_old = phi_new;
        
        % compute advection velocity at cell centers
        phi_x(:,1) = (phi_old(:,2) - phi_old(:,1))/dx;
        phi_x(:,2:end-1) = (phi_old(:,3:end) - phi_old(:,1:end-2))/(2*dx);
        phi_x(:,end) = (phi_old(:,end) - phi_old(:,end-1))/dx;
        phi_y(1,:) = (phi_old(2,:) - phi_old(1,:))/dx;
        phi_y(2:end-1,:) = (phi_old(3:end,:) - phi_old(1:end-2,:))/(2*dx);
        phi_y(end,:) = (phi_old(end,:) - phi_old(end-1,:))/dx;
        
        u = sgn_vals.*phi_x./sqrt(phi_x.^2 + phi_y.^2);
        v = sgn_vals.*phi_y./sqrt(phi_x.^2 + phi_y.^2);
        dtau = CFL*dx/max([max(abs(u)),max(abs(v))]); % adaptive timestep
        
        phi_new = advect2D(x,y,phi_old,u,v,dtau,sgn_vals);
        if (delay > 0)
            figure(1)
            plot3(x,y,psi,'ro',x,y,phi_new,'bs');
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