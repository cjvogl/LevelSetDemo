function psi_np1 = advect2D(x,y,psi_n,u_nph,v_nph,dt,f_nph)
%   solves the advection equation with extrapolation boundary conditions
%
%       psi_t(X,t) + u(X,t)*psi_x(X,t) + v(X,t)*psi_y(X,t) = f(X,t)
%
%   inputs:
%       x - meshgrid array containing x_j (assuming dx == dy)
%       y - meshgrid array containing y_i (assuming dy == dx)
%       psi_n - array containing psi(x_j,y_i,t_n)
%       u_nph - array containing either u(x_j,y_i,t_n + 0.5*dt) or
%               u & v at (x_{j\pm1/2},y_{i\pm1/2}, t_n + 0.5*dt)
%       dt - time step
%       f_nph - array containing f(x_j,y_i,t_n + 0.5*dt) (optional)
%
%   outputs:
%       psi_np1 - array containing psi(x_j,y_i,t_n + dt)

    if (nargin == 6)
        f_nph = zeros(size(x));
    end

    dx = x(1,2) - x(1,1);
  
    % add ghost nodes
    s = size(x);
    
    psi = zeros(s(1)+2,s(2)+2);
    psi(2:end-1,1) = 2.0*psi_n(:,1) - psi_n(:,2);
    psi(2:end-1,end) = 2.0*psi_n(:,end) - psi_n(:,end-1);
    psi(2:end-1,2:end-1) = psi_n;
    psi(end,2:end-1) = 2.0*psi_n(end,:) - psi_n(end-1,:);
    psi(1,2:end-1) = 2.0*psi_n(1,:) - psi_n(2,:);

    f = zeros(s(1)+2,s(2)+2);
    f(2:end-1,1) = 2.0*f_nph(:,1) - f_nph(:,2);
    f(2:end-1,end) = 2.0*f_nph(:,end) - f_nph(:,end-1);
    f(2:end-1,2:end-1) = f_nph;
    f(end,2:end-1) = 2.0*f_nph(end,:) - f_nph(end-1,:);
    f(1,2:end-1) = 2.0*f_nph(1,:) - f_nph(2,:);
    
    u = zeros(s(1)+2,s(1)+2);
    u(2:end-1,1) = 2.0*u_nph(:,1) - u_nph(:,2);
    u(2:end-1,end) = 2.0*u_nph(:,end) - u_nph(:,end-1);
    u(2:end-1,2:end-1) = u_nph;
    u(end,2:end-1) = 2.0*u_nph(end,:) - u_nph(end-1,:);
    u(1,2:end-1) = 2.0*u_nph(1,:) - u_nph(2,:);
    
    v = zeros(s(1)+2,s(1)+2);
    v(2:end-1,1) = 2.0*v_nph(:,1) - v_nph(:,2);
    v(2:end-1,end) = 2.0*v_nph(:,end) - v_nph(:,end-1);
    v(2:end-1,2:end-1) = v_nph;
    v(end,2:end-1) = 2.0*v_nph(end,:) - v_nph(end-1,:);
    v(1,2:end-1) = 2.0*v_nph(1,:) - v_nph(2,:);
    
    % compute derivatives
    dpsi_x_i_jmh = (psi(2:end-1,2:end-1) - psi(2:end-1,1:end-2))/dx;
    dpsi_x_i_jph = (psi(2:end-1,3:end) - psi(2:end-1,2:end-1))/dx;
    dpsi_y_imh_j = (psi(2:end-1,2:end-1) - psi(1:end-2,2:end-1))/dx;
    dpsi_y_iph_j = (psi(3:end,2:end-1) - psi(2:end-1,2:end-1))/dx;
    
    dpsi_xy_imh_jmh = (psi(2:end-1,2:end-1) - psi(2:end-1,1:end-2) ...
                       - psi(1:end-2,2:end-1) + psi(1:end-2,1:end-2))/dx^2;
    dpsi_xy_iph_jmh = (psi(3:end,2:end-1) - psi(3:end,1:end-2) ...
                       - psi(2:end-1,2:end-1) + psi(2:end-1,1:end-2))/dx^2;
    dpsi_xy_imh_jph = (psi(2:end-1,3:end) - psi(2:end-1,2:end-1) ...
                       - psi(1:end-2,3:end) + psi(1:end-2,2:end-1))/dx^2;
    dpsi_xy_iph_jph = (psi(3:end,3:end) - psi(3:end,2:end-1) ...
                       - psi(2:end-1,3:end) + psi(2:end-1,2:end-1))/dx^2;
    
    df_x_i_jmh = (f(2:end-1,2:end-1) - f(2:end-1,1:end-2))/dx;
    df_x_i_jph = (f(2:end-1,3:end) - f(2:end-1,2:end-1))/dx;
    df_y_imh_j = (f(2:end-1,2:end-1) - f(1:end-2,2:end-1))/dx;
    df_y_iph_j = (f(3:end,2:end-1) - f(2:end-1,2:end-1))/dx;

    u_imh_j = 0.5*(u(1:end-2,2:end-1) + u(2:end-1,2:end-1));
    u_iph_j = 0.5*(u(2:end-1,2:end-1) + u(3:end,2:end-1));
    u_i_jmh = 0.5*(u(2:end-1,1:end-2) + u(2:end-1,2:end-1));
    u_i_jph = 0.5*(u(2:end-1,2:end-1) + u(2:end-1,3:end));
    du_x_i_j = (u(2:end-1,3:end) - u(2:end-1,1:end-2))/(2*dx);
    du_y_i_j = (u(3:end,2:end-1) - u(1:end-2,2:end-1))/(2*dx);

    v_imh_j = 0.5*(v(1:end-2,2:end-1) + v(2:end-1,2:end-1));
    v_iph_j = 0.5*(v(2:end-1,2:end-1) + v(3:end,2:end-1));
    v_i_jmh = 0.5*(v(2:end-1,1:end-2) + v(2:end-1,2:end-1));
    v_i_jph = 0.5*(v(2:end-1,2:end-1) + v(2:end-1,3:end));
    dv_x_i_j = (v(2:end-1,3:end) - v(2:end-1,1:end-2))/(2*dx);
    dv_y_i_j = (v(3:end,2:end-1) - v(1:end-2,2:end-1))/(2*dx);

    
    % compute advection using scheme from C.J.Vogl, SISC, 2016
    u_plus = (sign(u_i_jmh) > 0).*u_i_jmh;
    u_minus = (sign(u_i_jph) < 0).*u_i_jph;
    v_plus = (sign(v_imh_j) > 0).*v_imh_j;
    v_minus = (sign(v_iph_j) < 0).*v_iph_j;
    psi_np1 = psi_n - dt*(1 + 0.5*dt*du_x_i_j).*u_plus.*dpsi_x_i_jmh ...
                    - dt*(1 + 0.5*dt*du_x_i_j).*u_minus.*dpsi_x_i_jph ...
                    - dt*(1 + 0.5*dt*dv_y_i_j).*v_plus.*dpsi_y_imh_j ...
                    - dt*(1 + 0.5*dt*dv_y_i_j).*v_minus.*dpsi_y_iph_j ...
        + 0.5*dt*abs(u_i_jmh).*(1 - dt/dx*abs(u_i_jmh)).*dpsi_x_i_jmh ...
        - 0.5*dt*abs(u_i_jph).*(1 - dt/dx*abs(u_i_jph)).*dpsi_x_i_jph ...
        + 0.5*dt*abs(v_imh_j).*(1 - dt/dx*abs(v_imh_j)).*dpsi_y_imh_j ...
        - 0.5*dt*abs(v_iph_j).*(1 - dt/dx*abs(v_iph_j)).*dpsi_y_iph_j ...
            + 0.5*dt^2*du_y_i_j.*v_i_jmh.*(u_plus > 0).*dpsi_x_i_jmh ...
            - 0.5*dt^2*du_y_i_j.*v_i_jph.*(u_minus < 0).*dpsi_x_i_jph ...
            + 0.5*dt^2*dv_x_i_j.*u_imh_j.*(v_plus > 0).*dpsi_y_imh_j ...
            - 0.5*dt^2*dv_x_i_j.*u_iph_j.*(v_minus < 0).*dpsi_y_iph_j ...
                    + dt^2*(u_plus.*v_plus.*dpsi_xy_imh_jmh ...
                            + u_plus.*v_minus.*dpsi_xy_iph_jmh ...
                            + u_minus.*v_plus.*dpsi_xy_imh_jph ...
                            + u_minus.*v_minus.*dpsi_xy_iph_jph) ...
                    + dt*f_nph ...
                    - 0.5*dt^2*(u_plus.*df_x_i_jmh + u_minus.*df_x_i_jph)...         
                    - 0.5*dt^2*(v_plus.*df_y_imh_j + v_minus.*df_y_iph_j);         
                
end
    
    