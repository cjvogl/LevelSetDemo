function phi = reinitializeFMM1D(x,psi,delay)
%   reinitializes the level set psi(x) into a signed-distance function:
%
%       phi_x(x,t)^2 = 1
%
%   using the fast marching method
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
    elseif (delay > 0)
        close all;
    end
    
    dx = x(2) - x(1);
    
    % initialize fast marching method variables
    accepted = zeros(size(x));
    trial = inf*ones(size(x)); % indicates distant node
    psi_interp = @(z)(interp1(x,psi,z,'pchip'));
    
    % mark initial accepted nodes
    recent_accepted = find(psi(2:end-1).*psi(3:end) <= 0 ...
                            |psi(2:end-1).*psi(1:end-2) <= 0)+1;
    if (psi(1)*psi(2) <= 0)
        recent_accepted(end+1) = 1;
    end
    if (psi(end)*psi(end-1) <= 0)
        recent_accepted(end+1) = length(x);
    end
    
    % seed values for accepted nodes
    for j=1:length(recent_accepted)
        current_ind = recent_accepted(j);
        x_star = fzero(psi_interp, x(current_ind));
        
        d = abs(x(current_ind) - x_star);
        accepted(current_ind) = d*sign(psi(current_ind));
        trial(current_ind) = -inf; % indicates accepted node
    end
    
    % complete fast marching reintialization
    while (max(trial) > -inf)
        if (delay > 0)
            figure(1)
            plot(x,accepted,'bo',x,psi,'--r');
            pause(delay);
        end
        
        % compute trial values for non-accepted neighbors of accepted nodes
        for j=1:length(recent_accepted)
            current_ind = recent_accepted(j);
            s = sign(psi(current_ind));
            if (current_ind > 1 && trial(current_ind-1) > -inf)
                trial(current_ind-1) = accepted(current_ind) + s*dx;
            end
            if (current_ind < length(x) && trial(current_ind+1) > -inf)
                trial(current_ind+1) = accepted(current_ind) + s*dx;
            end
        end
        
        % find trial node closest to interface to move into accepted
        recent_accepted = find(abs(trial) == min(abs(trial)),1);
        accepted(recent_accepted) = trial(recent_accepted);
        trial(recent_accepted) = -inf;
    end
                
    if (delay > 0)
        plot(x,accepted,'bo',x,psi,'--r');
        pause(delay);
    end
    
    phi = accepted;
    
end