% Driver file to simulate a moving interface, with position x(t) and 
% velocity u(x(t),t), in the domain [a,b] using the level set method
% (so that psi(x(t),t) = 0).

% user-supplied parameters
a = 0;                                  % left boundary of domain
b = 1;                                  % right boundary of domain
M = 11;                                % # of grid points
T = 5;                                 % final time
N = 50;                                % # of time steps to take
delay = 0.01;                           % animation delay (0: no animation)
u = @(x,t)(cos(2*pi*t)*ones(size(x)));  % domain velocity u(x,t)
phi0 = @(x)(atan((x-0.5*(a+b))*10));     % initial level set psi(x,0)

% set up domain and solution variables
x = linspace(a,b,M);
dt = T/N;
X = zeros(N+1,1);
t = 0:dt:T;

% create initial level set and save interface position
phi = phi0(x);
% phi = reinitializeFMM1D(x,phi,delay); %optional
% phi = reinitializePDE1D(x,phi,delay); %optional
X(1) = fzero(@(y)(interp1(x,phi,y,'pchip')),0.5*(a+b));

if (delay > 0)
    subplot(1,2,1)
    plot(x,phi,'-b',X(1),0,'ks','linewidth',2,'markersize',8)
    title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
    xlabel('x','fontsize',12,'fontweight','bold');
    ylabel('Level Set','fontsize',12,'fontweight','bold');
    subplot(1,2,2)
    plot(t(1:1),X(1:1),'-ok','linewidth',2,'markersize',8)
    title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
    xlabel('t','fontsize',12,'fontweight','bold');
    ylabel('Interface Position','fontsize',12,'fontweight','bold');
    pause;
end

% main solution loop
for j=1:N
    
    % compute velocity at t_n + 0.5*dt
    u_nph = u(x,t(j)+0.5*dt);
    
    % advect level set
    phi = advect1D(x,phi,u_nph,dt);
    
    % reinitialize level set (optional)
    % phi = reinitializeFMM1D(x,phi);
    % phi = reinitializePDE1D(x,phi);
    
    % save interface position
    X(j+1) = fzero(@(y)(interp1(x,phi,y)),X(j),'pchip');
   
    if (delay > 0)
        subplot(1,2,1)
        plot(x,phi,'-b',X(j+1),0,'ks','linewidth',2,'markersize',8)
        title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
        xlabel('x','fontsize',12,'fontweight','bold');
        ylabel('Level Set','fontsize',12,'fontweight','bold');
        subplot(1,2,2)
        plot(t(1:j+1),X(1:j+1),'-ok','linewidth',2,'markersize',8)
        title(sprintf('time %f', t),'fontsize',12,'fontweight','bold');
        xlabel('t','fontsize',12,'fontweight','bold');
        ylabel('Interface Position','fontsize',12,'fontweight','bold');
        pause(delay);
    end
    
end