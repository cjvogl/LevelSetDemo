function phi = reinitializeFMM2D(x,y,psi,delay)
%   reinitializes the level set psi(x) into a signed-distance function:
%
%       phi_x(x,y,t)^2 + phi_y(x,y,t)^2 = 1
%
%   using the fast marching method
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
    end
    
    dx = x(1,2) - x(1,1);
    
    % initialize fast marching method variables
    accepted = zeros(size(x));
    trial = inf*ones(size(x)); % indicates distant node
    psi_interp = @(X)(interp2(x,y,psi,X(1),X(2),'cubic'));
    
    % mark initial accepted nodes
    tmpInd = find(psi(2:end-1,2:end-1).*psi(3:end,2:end-1) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(1:end-2,2:end-1) <= 0 ... 
                    | psi(2:end-1,2:end-1).*psi(2:end-1,3:end) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(2:end-1,1:end-2) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(3:end,3:end) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(1:end-2,3:end) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(3:end,1:end-2) <= 0 ...
                    | psi(2:end-1,2:end-1).*psi(1:end-2,1:end-2) <= 0);
    [tmpI,tmpJ] = ind2sub([length(x)-2,length(x)-2],tmpInd);
    recent_accepted = [tmpI+1,tmpJ+1];
    if (psi(1,1)*psi(1,2) <= 0 || psi(1,1)*psi(2,2) <= 0 ...
            || psi(1,1)*psi(2,1) <= 0)
        recent_accepted(end+1,:) = [1,1];
    end
    if (psi(1,end)*psi(1,end-1) <= 0 || psi(1,end)*psi(2,end-1) <= 0 ...
            || psi(1,end)*psi(2,end) <= 0)
        recent_accepted(end+1,:) = [1,length(x)];
    end
    for i=2:length(x)-1
        if (psi(i,1)*psi(i-1,1) <= 0 || psi(i,1)*psi(i-1,2) <= 0 ...
                || psi(i,1)*psi(i,2) <= 0 || psi(i,1)*psi(i+1,2) <= 0 ...
                || psi(i,1)*psi(i+1,1) <= 0)
            recent_accepted(end+1,:) = [i,1];
        end
        if (psi(i,end)*psi(i-1,end) <= 0 || psi(i,end)*psi(i-1,end-1) <= 0 ...
                || psi(i,end)*psi(i,end-1) <= 0 || psi(i,end)*psi(i+1,end-1) <= 0 ...
                || psi(i,end)*psi(i+1,end) <= 0)
            recent_accepted(end+1,:) = [i,length(x)];
        end
    end
    if (psi(end,1)*psi(end-1,1) <= 0 || psi(end,1)*psi(end-1,2) <= 0 ...
            || psi(end,1)*psi(end,2) <= 0)
        recent_accepted(end+1,:) = [length(x),1];
    end
    if (psi(end,end)*psi(end-1,end) <= 0 || psi(end,end)*psi(end-1,end-1) <= 0 ...
            || psi(end,end)*psi(end,end-1) <= 0)
        recent_accepted(end+1,:) = [length(x),length(x)];
    end
    for j=2:length(x)-1
        if (psi(1,j)*psi(1,j-1) <= 0 || psi(1,j)*psi(2,j-1) <= 0 ...
                || psi(1,j)*psi(2,j) <= 0 || psi(1,j)*psi(2,j+1) <= 0 ...
                || psi(1,j)*psi(1,j+1) <= 0)
            recent_accepted(end+1,:) = [1,j];
        end
        if (psi(end,j)*psi(end,j-1) <= 0 || psi(end,j)*psi(end-1,j-1) <= 0 ...
                || psi(end,j)*psi(end-1,j) <= 0 || psi(end,j)*psi(end-1,j+1) <= 0 ...
                || psi(end,j)*psi(end,j+1) <= 0)
            recent_accepted(end+1,:) = [length(x),j];
        end
    end
    
    % seed values for accepted nodes
    for j=1:size(recent_accepted,1)
        currentI = recent_accepted(j,1);
        currentJ = recent_accepted(j,2);
        currentX = x(currentI,currentJ);
        currentY = y(currentI,currentJ);
        
        % NOTE that a constrained optimization routine is usually used
        % here, but in the interest of keeping this code accessible, the
        % following langrange multiplier/penalty is used
        lambda = 1e6;
        obj_func = @(Xs)((currentX-Xs(1)).^2 + (currentY-Xs(2)).^2 + lambda*(psi_interp(Xs))^2);
        X_star = fminsearch(obj_func,[currentX,currentY]);
        
        % If you have the optimization toolbox, it is highly recommended
        % that you use fmincon instead
        %obj_func = @(Xs)((currentX-Xs(1)).^2 + (currentY-Xs(2)).^2);
        %con_func = @(X)(level_set_constraint(X,psi_interp));
        %options = optimoptions('fmincon','Algorith','active-set','Display','off');
        %X_star = fmincon(obj_func,[currentX,currentY],[],[],[],[],...
        %            [x(1,1),y(1,1)],[x(end,end),y(end,end)],con_func,options);

        x_star = X_star(1);
        y_star = X_star(2);
        
        d = sqrt((currentX - x_star)^2 + (currentY - y_star)^2);
        accepted(currentI,currentJ) = d*sign(psi(currentI,currentJ));
        trial(currentI,currentJ) = -inf; % indicates accepted node
    end
    
    % complete fast marching reintialization
    while (max(max(trial)) > -inf)
        if (delay > 0)
            figure(1)
            plot3(x,y,psi,'ro',x,y,accepted,'bs');
            pause(delay);
        end
        
        % compute trial values for non-accepted neighbors of accepted nodes
        for j=1:size(recent_accepted,1)
            currentI = recent_accepted(j,1);
            currentJ = recent_accepted(j,2);
            s = sign(psi(currentI,currentJ));
            if (currentI > 1 && trial(currentI-1,currentJ) > -inf)
                trial(currentI-1,currentJ) = obtainTrialValue(currentI-1,...
                                    currentJ,accepted,trial,dx,s);
            end
            if (currentI < length(x) && trial(currentI+1,currentJ) > -inf)
                trial(currentI+1,currentJ) = obtainTrialValue(currentI+1,...
                                    currentJ,accepted,trial,dx,s);
            end
            if (currentJ > 1 && trial(currentI,currentJ-1) > -inf)
                trial(currentI,currentJ-1) = obtainTrialValue(currentI,...
                                    currentJ-1,accepted,trial,dx,s);
            end
            if (currentJ < length(x) && trial(currentI,currentJ+1) > -inf)
                trial(currentI,currentJ+1) = obtainTrialValue(currentI,...
                                    currentJ+1,accepted,trial,dx,s);
            end
        end
        
        % find trial node closest to interface to move into accepted
        tmpInd = find(abs(trial) == min(min(abs(trial))),1);
        accepted(tmpInd) = trial(tmpInd);
        trial(tmpInd) = -inf;
        [tmpI,tmpJ] = ind2sub(size(x),tmpInd);
        recent_accepted = [tmpI,tmpJ];
    end
                
    if (delay > 0)
        figure(1)
        plot3(x,y,psi,'ro',x,y,accepted,'bs');
        pause(delay);
    end
    
    phi = accepted;
    
end

function phi = obtainTrialValue(i,j,accepted,trial,dx,s)
%   function to generate a trial signed-distance level set value for 
%   (x(i,j),y(i,j)) by extending the information from accepted nodes.
%   Upwind directions are promoted so that information propagates in the 
%   correct direction

    % First, extension from 2 neighboring accepted nodes is attempted in
    % the x direction, using y_{i-1}
    if (i > 1)
        phi_b = attemptLeftAndRightExtension(i,j,-1,trial,accepted,dx,s);
        dphi_y_b = (phi_b - accepted(i-1,j))/dx;
    else
        phi_b = inf;
    end
    
    % Now, extension from 2 neighboring accepted nodes is attempted in
    % the x direction, using y_{i+1}
    if (i < length(accepted))
        phi_t = attemptLeftAndRightExtension(i,j,1,trial,accepted,dx,s);
        dphi_y_t = (accepted(i+1,j) - phi_t)/dx;
    else
        phi_t = inf;
    end

    % Promote the more upwind of the valid extensions
    if (phi_b < inf && phi_t < inf)
        if (abs(dphi_y_b) > abs(dphi_y_t))
            phi = phi_b;
        else
            phi = phi_t;
        end
    elseif (phi_b < inf)
        phi = phi_b;
    elseif (phi_t < inf)
        phi = phi_t;
    else
        phi = inf;
    end    
    
    % If all 2 point extensions are invalid, use 1 point extension
    if (phi == inf)
        phi_l = inf;
        phi_r = inf;
        phi_b = inf;
        phi_t = inf;
        if (j > 1 && trial(i,j-1) == -inf)
            phi_l = accepted(i,j-1) + s*dx;
        end
        if (j < length(accepted) && trial(i,j+1) == -inf)
            phi_r = accepted(i,j+1) + s*dx;
        end
        if (i > 1 && trial(i-1,j) == -inf)
            phi_b = accepted(i-1,j) + s*dx;
        end
        if (i < length(accepted) && trial(i+1,j) == -inf)
            phi_t = accepted(i+1,j) + s*dx;
        end
        phi = s*min(abs([phi_l,phi_r,phi_b,phi_t]));
    end

end

function phi = attemptLeftAndRightExtension(i,j,d,trial,accepted,dx,s)
%   function to generate a trial signed-distance level set value for 
%   (x(i,j),y(i,j)) by extending the information from accepted nodes to the
%   left and right, using the bottom neighbor if d = -1 or the top
%   neighbor if d = 1.  Upwind directions are promoted.

    % initialize trial values to infinity
    phi_l = inf;
    phi_r = inf;
    dphi_x_l = 0;
    dphi_x_r = 0;
    
    % Obtain all valid left and right extensions
    if (j > 1 && trial(i+d,j) == -inf && trial(i,j-1) == -inf)
        A = accepted(i+d,j);
        B = accepted(i,j-1);
        if (2 - ((A-B)/dx)^2 >= 0)
            phi_l = 0.5*(A+B) + 0.5*dx*s*sqrt(2 - ((A-B)/dx)^2);
            dphi_x_l = (phi_l - B)/dx;
        end
        if (phi_l*s < A*s || phi_l*s < B*s) % check for invalid extension
            phi_l = inf;
        end
    end
    if (j < length(accepted) && trial(i+d,j) == -inf && trial(i,j+1) == -inf)
        A = accepted(i+d,j);
        B = accepted(i,j+1);
        if (2 - ((A-B)/dx)^2 >= 0)
            phi_r = 0.5*(A+B) + 0.5*dx*s*sqrt(2 - ((A-B)/dx)^2);
            dphi_x_r = (B - phi_r)/dx;
        end
        if (phi_r*s < A*s || phi_r*s < B*s) % check for invalid extension
            phi_r = inf;
        end
    end
    
    % Promote the more upwind of the valid extensions
    if (phi_l < inf && phi_r < inf)
        if (abs(dphi_x_l) > abs(dphi_x_r))
            phi = phi_l;
        else
            phi = phi_r;
        end
    elseif (phi_l < inf)
        phi = phi_l;
    elseif (phi_r < inf)
        phi = phi_r;
    else
        phi = inf;
    end
end

function [c,ceq] = level_set_constraint(X,psi_interp)
%   function to be used in constrained miminimzation if available
%   requires that psi(X(1),X(2)) = 0 for cubic interpolant of psi.

    c = [];
    ceq = psi_interp(X);
end