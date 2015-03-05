function Project2c(alpha, n_pan, m, p, tt)
    % Project2jfd - 2D potential flow via source panel method
    %    written by John Dannenhoffer for AEE342 (Spring 2015)
    %    adapted from Anderson's "Fundamentals of Aerodynamics, 5th ed"
    %    edited by Joel Lubinitsky for airfoil
    %
    % alpha  = angle of attack        (deg)
    % series = NACA series designator (mptt)
    % n_pan  = number of panels

    close all
    
    n_sl  =   5;                % number of streamlines
    dtbub = 0.25;                % time increment between bubbles
    xmin  = -0.5;                % mininum X in domain
    xmax  = +1.5;                % maximum X in domain
    ymin  = -3.0;                % minimum Y in domain
    ymax  = +3.0;                % maximum Y in domain

    % format inputs
    n_pan = n_pan + 2;
    m = m * 0.01;
    p = p * 0.1;
    tt = tt * 0.01;

    % define airfoil geometry
    zeta = linspace(pi, 2 * pi, n_pan / 2);
    xCamber = 0.5 * (1 + cos(zeta));
%     yThickness = (tt / 0.20) * ((0.2969 * sqrt(xCamber)) - (0.1260 .* xCamber) - (0.3516 .* xCamber .^ 2) + (0.2843 .* xCamber .^ 3) - (0.1015 .* xCamber .^ 4));
    yThickness = (tt / 0.20) * ((0.2969 * sqrt(xCamber)) - (0.1260 .* xCamber) - (0.3516 .* xCamber .^ 2) + (0.2843 .* xCamber .^ 3) - (0.1036 .* xCamber .^ 4));

    if p > 0 && m > 0
    indexP = find(xCamber < p, 1, 'last');

    xCamber1 = xCamber(1 : indexP);
    xCamber2 = xCamber(indexP + 1 : end);
    yCamber1 = (m / p ^ 2) * (2 * p * xCamber1 - xCamber1 .^ 2);
    yCamber2 = (m / (1 - p) ^ 2) * (1 - (2 * p) + (2 * p * xCamber2) - xCamber2 .^ 2);
    yCamber = [yCamber1, yCamber2];

    dydxCamber1 = (2 * m / p) * (1 - xCamber1 / p);
    dydxCamber2 = (2 * m / (1 - p) ^ 2) * (p - xCamber2);
    dydxCamber = [dydxCamber1, dydxCamber2];
    theta = atan(dydxCamber);

    xUpper = xCamber - yThickness .* sin(theta);
    xUpper = xUpper(1 : 1 : end - 1);
    xLower = xCamber + yThickness .* sin(theta);
    xLower = xLower(end : -1 : 2);

    yUpper = yCamber + yThickness .* cos(theta);
    yUpper = yUpper(1 : 1 : end - 1);
    yLower = yCamber - yThickness .* cos(theta);
    yLower = yLower(end : -1 : 2);

    X = [xUpper, xLower];
    Y = [yUpper, yLower];
    
    elseif p == 0 && m == 0
        xUpper = xCamber;
        xLower = xCamber;
        
        yUpper = yThickness;
        yLower = -yThickness;
        
        X = [xUpper(1 : 1 : end - 1), xLower(end : -1 : 2)];
        Y = [yUpper(1 : 1 : end - 1), yLower(end : -1 : 2)];
        
    end

    % plot the configuration with first point repeated (figure 1)
    figure(1)
    hold on
    plot([X, X(1)], [Y, Y(1)], '-*')
%     plot(X(500), Y(500), 'o', 'color', [1 0 0])
    title(strcat('Original configuration - NACA ', num2str(m * 100), num2str(p * 10), num2str(tt * 100)))
    xlabel('x')
    ylabel('y')
    axis equal
    grid on

    % define the freestream velocities in the x- and y- directions
    Uinf  = cos(alpha * pi/180);
    Vinf  = sin(alpha * pi/180);
    


    % find the source panel strengths
    [xcp, ycp, Cp, lambda] = sourcePanel(Uinf, Vinf, X, Y);
    
        % Superposition of vortex at LE
    strengthGamma = 0;
    
%     tolerKutta = 0.01;
%     step = 0.01;
%     
%     while sqrt(Uvel(1, 0, Uinf, X, Y, lambda, strengthGamma) ^ 2 + Vvel(1, 0, Uinf, X, Y, lambda, strengthGamma) ^ 2) > tolerKutta
%         if Vvel(1, 0, Uinf, X, Y, lambda, strengthGamma) < 0
%             
    
    % integrate pressure coefficient for force coefficients
    CpUpper = Cp(1 : round(n_pan / 2) - 1);
    CpLower = Cp(round(n_pan / 2) : end);
    CpLower = CpLower(end : -1 : 1);
    
    xCpUpper = xcp(1 : round(n_pan / 2) - 1);
    xCpLower = xcp(round(n_pan / 2) : end);
    xCpLower = xCpLower(end : -1 : 1);
    
    Cn = Project2aCn(n_pan, xCpUpper, xCpLower, CpUpper, CpLower);
    Ca = Project2aCa(n_pan, xCpUpper, xCpLower, CpUpper, CpLower, yUpper, yLower);
    
    Cl = Cn * cos(alpha) - Ca * sin(alpha);
    Cd = Cn * sin(alpha) + Ca * cos(alpha);

    % plot the source panel strengths (figure 2, subplot 1)
    figure(2)
    subplot(2,2,1)
        plot(1:length(lambda), lambda, '-o')
    
        title(strcat('AEE342 - Project2a, \alpha = ', num2str(alpha), '^\circ'))
        xlabel('Panel number')
        ylabel('Source panel strength (lambda)')
        axis([0 140 -20 20])
        grid on
    
    % plots the surface Cp distribution (figure 2, subplots 2 and 3)
    subplot(2,2,2)
        plot([xcp, xcp(1)], -[Cp, Cp(1)])
    
        title(strcat('C_n =', num2str(Cn)))
        xlabel('x (control points)')
        ylabel('-Cp')
        xlim([-0.5 1.5])
        ylim([-1 4])
        grid on
    
    subplot(2,2,3)
        plot(-[Cp, Cp(1)], [ycp, ycp(1)])
    
        title(strcat('C_a =', num2str(Ca)))
        xlabel('-Cp')
        ylabel('y (control points)')
        ylim([-0.8 0.8])
        xlim([-1 4])
        grid on

    % plot of configuration and streamlines (figure 2, subplot 4)
    subplot(2,2,4)
        plot([X, X(1)], [Y, Y(1)], '-*b')
    
        hold on
        for i_sl = 1 : n_sl
            xbeg = -0.5;
            ybeg = -0.8 + (i_sl-1) / (n_sl-1) * (0.8 - -0.8);
        
            [t_sl, xy_sl] = ode45(@(t, xy) strmLine(t, xy, Uinf, Vinf, ...
                                    X, Y, lambda, strengthGamma), [0, 10], [xbeg, ybeg]);
            plot(xy_sl(:,1), xy_sl(:,2), '-g')

            % add bubbles to the streamlines (if desired)
            if (dtbub > 0)
                t_bub = 0 : dtbub : t_sl(end);
                x_bub = interp1(t_sl, xy_sl(:,1), t_bub);
                y_bub = interp1(t_sl, xy_sl(:,2), t_bub);
                plot(x_bub, y_bub, 'og', 'MarkerSize', 2)
            end % if
        end % for i_sl
        hold off
    
        title(strcat('C_l =', num2str(Cl), ' C_d =', num2str(Cd)))
        xlabel('x')
        ylabel('y')
        axis([-0.5 1.5 -0.8 0.8])
        
        

    
end % function Project2jfd

%--------------------------------------------------------------------

function [xcp, ycp, Cp, lambda] = sourcePanel(Uinf, Vinf, X, Y)
    % sourcePanel - source panel method in two dimensions
    %    written by John Dannenhoffer
    %
    % inputs:
    %    Uinf       x-component of freestream velocity
    %    Vinf       y-component of freestream velocity
    %    X          x-coordinate of boundary points
    %    Y          y-coordinate of boundary points
    % outputs:
    %    xcp        x-coordinate of control points
    %    ycp        y-coordinate of control points
    %    Cp         pressure coefficient at control points
    %    lambda     column-vector of panel strengths

    toler = 0.000000000000000000000001;                      % minimum panel length

    n = length(X);

    % find the control point locations (xcp, ycp), panel lengths (s),
    %    and panel inclinations (phi)
    xcp = zeros(1,n);
    ycp = zeros(1,n);
    S   = zeros(1,n);
    phi = zeros(1,n);

    for i = 1 : n
        if (i < n)
            ip1 = i + 1;
        else
            ip1 = 1;
        end % if
    
        xcp(i) =                           (X(ip1) + X(i))/2;
        ycp(i) =       (Y(ip1) + Y(i))/2;
        S(  i) = sqrt( (Y(ip1) - Y(i))^2 + (X(ip1) - X(i))^2);
        phi(i) = atan2((Y(ip1) - Y(i))   , (X(ip1) - X(i))  );
    
        % check that the points are distinct
        if (S(i) < toler)
            error('Points %d and %d are not distinct', i, ip1)
        end % if
    end % for i

    % set up the Mnorm matrix (which is called "I" in Anderson) and
    %    Mtang matrix (which is called "I'" in Anderson)
    Mnorm = zeros(n + 1, n + 1);
    Mtang = zeros(n, n);
    RHS   = zeros(n + 1, 1);

    % normal velocity at ith control point
    for i = 1 : n
        RHS(i) = -2 * pi * (Vinf * cos(phi(i)) - Uinf * sin(phi(i)));
    
        % contribution of jth panel
        for j = 1 : n
            if (i == j)
                Mnorm(i,j) = pi;
                Mtang(i,j) = 0;
            else
                A = - (xcp(i) - X(j)) * cos(phi(j)) - (ycp(i) - Y(j)) * sin(phi(j));
                B =   (xcp(i) - X(j)) ^ 2           + (ycp(i) - Y(j)) ^ 2;
                D = - (xcp(i) - X(j)) * sin(phi(i)) + (ycp(i) - Y(j)) * cos(phi(i));
                E =   (xcp(i) - X(j)) * sin(phi(j)) - (ycp(i) - Y(j)) * cos(phi(j));
                C = sin(phi(i) - phi(j));
        
                Mnorm(i,j) =    C/2        * log((S(j)^2 + 2*A*S(j) + B) / B) ...
                           + (D-A*C)/   E  * (atan((S(j)+A)/E) - atan(A/E));
                Mtang(i,j) = (D-A*C)/(2*E) * log((S(j)^2 + 2*A*S(j) + B) / B) ...
                           -    C          * (atan((S(j)+A)/E) - atan(A/E));
            end % if
        end % for j
    end % for i

    % solve for the source strengths
    lambda = Mnorm \ RHS;

    % compute the tangential velocities and hence the Cp
    V = Uinf*cos(phi) + Vinf*sin(phi);
    for j = 1 : n
        V = V + lambda(j) ./ (2*pi) .* Mtang(:,j)';
    end % for
    Cp = 1 - V.^2;
    
    figure(3)
    hold on
    axis([0 1 -1.5 1.5])
    plot(X, V)
    plot([0: 0.01 : 1], 0, 'color', [1 0 0])
    indexStagnation = find(V >= 0, 1, 'last')
    xStagnation = X(indexStagnation)
    plot(X(indexStagnation), V(indexStagnation), 'o', 'color', [1 0 0])

    check = sum(S*lambda);
    if (abs(check) > 1.0e-10)
        fprintf(1, 'sum of sources = %f (should be 0)\n', check);
    end % if
end % function sourcePanel

%--------------------------------------------------------------------

function u = Uvel(x, y, Uinf, X, Y, lambda, strengthGamma)
    % Uvel - x-component of velocity at (x,y)

    % uniform flow
    u = Uinf * ones(size(x));

    % loop through source panels
    n = length(lambda);
    for j = 1 : n
        if (j < n)
            jp1 = j + 1;
        else
            jp1 = 1;
        end % if
    
        S   = sqrt( (Y(jp1) - Y(j)).^2 + (X(jp1) - X(j)).^2);
        phi = atan2((Y(jp1) - Y(j))    , (X(jp1) - X(j))   );

        A = - (x - X(j)) .* cos(phi)  - (y - Y(j)) .* sin(phi);
        B =   (x - X(j)) .^ 2         + (y - Y(j)) .^ 2;
        D =   (x - X(j));
        E =   (x - X(j)) .* sin(phi)  - (y - Y(j)) .* cos(phi);
        C = sin(-pi/2 - phi);
        
        u = u + lambda(j)/(2*pi) * (    C/2     .* log((S^2 + 2*A*S + B) ./ B) ...
                                   + (D-A*C)./E .* (atan((S+A)./E) - atan(A./E)))...
                                   + (strengthGamma ./ (2 .* pi .* sqrt(x .^ 2 + y .^ 2)) .* sin(atan2(y, -x)));
    end % for j
end % function Uvel

%--------------------------------------------------------------------

function v = Vvel(x, y, Vinf, X, Y, lambda, strengthGamma)
    % Vvel - y-component of velocity at (x,y)

    % uniform flow
    v = Vinf * ones(size(x));

    % loop through source panels
    n = length(lambda);
    for j = 1 : n
        if (j < n)
            jp1 = j + 1;
        else
            jp1 = 1;
        end % if
    
        S   = sqrt( (Y(jp1) - Y(j)).^2 + (X(jp1) - X(j)).^2);
        phi = atan2((Y(jp1) - Y(j))    , (X(jp1) - X(j))  );

        A = - (x - X(j)) .* cos(phi)  - (y - Y(j)) .* sin(phi);
        B =   (x - X(j)) .^ 2         + (y - Y(j)) .^ 2;
        D =                           + (y - Y(j));
        E =   (x - X(j)) .* sin(phi)  - (y - Y(j)) .* cos(phi);
        C = sin(- phi);
        
        v = v + lambda(j)/(2*pi) * (    C/2     .* log((S^2 + 2*A*S + B) ./ B) ...
                                   + (D-A*C)./E .* (atan((S+A)./E) - atan(A./E)))...
                                   + (strengthGamma ./ (2 .* pi .* sqrt(x .^ 2 + y .^ 2)) .* cos(atan2(y, -x)));
    end % for j
end % function Vvel

%--------------------------------------------------------------------

function uv = strmLine(t, xy, Uinf, Vinf, X, Y, lambda, strengthGamma)
    % strmLine - vector velocity function (used in ode45)

    % extract the x and y coordinates
    x = xy(1);
    y = xy(2);

    % get the velocities and return them as a vector
    uv = [Uvel(x, y, Uinf, X, Y, lambda, strengthGamma); ...
         Vvel(x, y, Vinf, X, Y, lambda, strengthGamma)];
end % function strmLine
  
%--------------------------------------------------------------------
