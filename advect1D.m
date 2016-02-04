function psi_np1 = advect1D(x,psi_n,u_nph,dt,f_nph)
%   solves the advection equation with extrapolation boundary conditions
%
%       psi_t(x,t) + u(x,t)*psi_x(x,t) = f(x,t)
%
%   inputs:
%       x - vector containing x_i
%       psi_n - vector containing psi(x_i,t_n)
%       u_nph - vector containing either u(x_i,t_n + 0.5*dt) or
%               u(x_{i\pm1/2}, t_n + 0.5*dt)
%       dt - time step
%       f_nph - vector containing f(x_i,t_n + 0.5*dt) (optional)
%
%   outputs:
%       psi_np1 - vector containing psi(x_i,t_n + dt)

    if (nargin == 4)
        f_nph = zeros(size(x));
    end

    dx = x(2) - x(1);

    % add ghost nodes for psi and f
    s = size(x);
    s(s > 1) = s(s > 1) + 2;
    
    psi = zeros(s);
    psi(1) = 2.0*psi_n(1) - psi_n(2);
    psi(2:end-1) = psi_n;
    psi(end) = 2.0*psi_n(end) - psi_n(end-1);

    f = zeros(s);
    f(1) = 2.0*f_nph(1) - f_nph(2);
    f(2:end-1) = f_nph;
    f(end) = 2.0*f_nph(end) - f_nph(end-1);    
    
    % set edge velocities and derivatives (computing averages if necesary)
    if (length(u_nph) > length(x))
        u_imh = u_nph(1:end-1);
        u_iph = u_nph(2:end);
        du_i = (u_nph(2:end) - u_nph(1:end-1))/dx;
    else
        u = zeros(s); % add ghost nodes first
        u(1) = 2.0*u_nph(1) - u_nph(2);
        u(2:end-1) = u_nph;
        u(end) = 2.0*u_nph(end) - u_nph(end-1);
        u_imh = 0.5*(u(1:end-2) + u(2:end-1));
        u_iph = 0.5*(u(2:end-1) + u(3:end));
        du_i = (u(3:end) - u(1:end-2))/(2*dx);
    end
    
    % compute derivatives for psi and f
    dpsi_imh = (psi(2:end-1) - psi(1:end-2))/dx;
    dpsi_iph = (psi(3:end) - psi(2:end-1))/dx;
    df_imh = (f(2:end-1) - f(1:end-2))/dx;
    df_iph = (f(3:end) - f(2:end-1))/dx;

    
    % compute advection using scheme from C.J.Vogl, SISC, 2016
    u_plus = (sign(u_imh) > 0).*u_imh;
    u_minus = (sign(u_iph) < 0).*u_iph;
    psi_np1 = psi_n - dt*(1 + 0.5*dt*du_i).*u_plus.*dpsi_imh ...
                    - dt*(1 + 0.5*dt*du_i).*u_minus.*dpsi_iph ...
                    + 0.5*dt*abs(u_imh).*(1 - dt/dx*abs(u_imh)).*dpsi_imh ...
                    - 0.5*dt*abs(u_iph).*(1 - dt/dx*abs(u_iph)).*dpsi_iph ...
                    + dt*f_nph ...
                    - 0.5*dt^2*(u_plus.*df_imh + u_minus.*df_iph);         
                
end
    
    