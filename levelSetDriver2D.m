% Driver file to simulate a moving interface, parameterized by 
% (x(s,t),y(s,t)) with velocity [u(x(s,t),y(s,t),t),v(x(s,t),y(s,t),t)],
% in the domain [a,b]x[a,b] using the level set method 
% (so that psi(x(s,t),y(s,t),t) = 0).

% user-supplied parameters
a = 0;                                  % left/bottom boundary of domain
b = 1;                                  % right/top boundary of domain
M = 21;                                % # of grid points
T = 5;                                 % final time
N = 100;                                % # of time steps to take
delay = 0.01;                           % animation delay (0: no animation)
u = @(x,y,t)(cos(2*pi*t)*ones(size(x)));  % horiz. domain velocity u(x,y,t)
v = @(x,y,t)(zeros(size(x)));             % vert. domain velocity v(x,y,t)
phi0 = @(x,y)((x-0.5*(a+b)).^2 + (y-0.5*(a+b)).^2 ...
                    - (0.15*(b-a))^2);     % initial level set psi(x,0)

% set up domain and solution variables
[x,y] = meshgrid(linspace(a,b,M),linspace(a,b,M));
dx = x(1,2)-x(1,1);
dt = T/N;
Phi = zeros(size(x,1),size(x,2),N+1);
t = 0:dt:T;

% create initial level set and velocity
phi = phi0(x,y);
% phi = reinitializeFMM2D(x,y,phi,delay); %optional
 phi = reinitializePDE2D(x,y,phi,delay); %optional
Phi(:,:,1) = phi;

% find max velocity to scale quiver plots
max_vel = 0;
for j=1:N
    u_temp = u(x,y,(j-0.5)*dt);
    v_temp = v(x,y,(j-0.5)*dt);
    max_vel = max([max_vel ...
                    max(max(sqrt(u_temp.*u_temp + v_temp.*v_temp)))]);
end
u_temp = u(x,y,0);
v_temp = v(x,y,0);
scale = dx/max_vel;

if (delay > 0)
    subplot(1,2,1)
    surfc(x,y,phi)
    title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
    xlabel('x','fontsize',12,'fontweight','bold');
    ylabel('y','fontsize',12,'fontweight','bold');
    zlabel('Level Set','fontsize',12,'fontweight','bold');
    subplot(1,2,2)
    contour(x,y,phi,[0 0],'linewidth',2)
    hold on;
    quiver(x,y,scale*u_temp,scale*v_temp,0);
    hold off;
    title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
    xlabel('x','fontsize',12,'fontweight','bold');
    ylabel('y','fontsize',12,'fontweight','bold');
    pause;
end

% main solution loop
for j=1:N
    
    % compute velocity at t_n + 0.5*dt
    u_nph = u(x,y,t(j)+0.5*dt);
    v_nph = v(x,y,t(j)+0.5*dt);
    
    % advect level set
    phi = advect2D(x,y,phi,u_nph,v_nph,dt);
    
    % reinitialize level set (optional)
    % phi = reinitializeFMM2D(x,phi);
    % phi = reinitializePDE2D(x,phi);
    
    % save interface position
    Phi(:,:,j+1) = phi;
   
    if (delay > 0)
        subplot(1,2,1)
        surfc(x,y,phi)
        title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
        xlabel('x','fontsize',12,'fontweight','bold');
        ylabel('y','fontsize',12,'fontweight','bold');
        zlabel('Level Set','fontsize',12,'fontweight','bold');
        subplot(1,2,2)
        contour(x,y,phi,[0 0],'linewidth',2)
        hold on;
        quiver(x,y,scale*u_nph,scale*v_nph,0)
        hold off;
        title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
        xlabel('x','fontsize',12,'fontweight','bold');
        ylabel('y','fontsize',12,'fontweight','bold');
        pause(delay);
    end
    
end