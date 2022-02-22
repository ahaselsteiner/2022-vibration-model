SYSTEM = 'A';

DO_ANIMATE = 0;
DO_WRITE_GIF = 0;

% System A, the table top experiment.
A.m1 = 0.0093; % in kg
A.m2 = 0.044; % in kg
A.k1 = 4.5; % in N/m
A.k2 = 0.001 * 180 / pi; % in Nm/rad
A.d = 0.038; % in m
A.I = 0.00005; % in kg m^2
A.x0 = 0.1;
A.y0 = 0.1;
A.xdot0 = 0;
A.ydot0 = 0;
A.max_time = 20;
A.label = 'tabletop';


% System B, the typical offshore wind turbine.
B.m1 = 800 * 10^3;
B.m2 = 400 * 10^3;
B.k1 = 5 * 10^6;
B.k2 = 3 * 10^9;
B.d = 0.3;
B.I = 4 * 10^7;
B.x0 = 1;
B.y0 = 1;
B.xdot0 = 1;
B.ydot0 = 0;
B.max_time = 120;
B.label = 'turbine_typical';

% System C, the unfavorable offshore wind turbine.
FACTOR_TORSION = 15;
FACTOR_DISTANCE = 3.33333;
FACTOR_IZZ = 1;
C.m1 = B.m1;
C.m2 = B.m2;
C.k1 = B.k1;
C.k2 = 1 / FACTOR_TORSION * B.k2;
C.d = FACTOR_DISTANCE * B.d;
C.I = FACTOR_IZZ * B.I;
C.x0 = B.x0;
C.y0 = B.y0;
C.xdot0 = B.xdot0;
C.ydot0 = B.ydot0;
C.max_time = B.max_time;
C.label = 'turbine_unfavorable';


damp_theta = 0.0001; % 0.0001 for  numerical stability (problems with system A otherwise).

switch SYSTEM
    case 'A'
        S = A;
    case 'B'
        S = B;
    case 'C'
        S = C;
    otherwise
        warning('Unexpected system.')
end
m1 = S.m1;
m2 = S.m2;
k1 = S.k1;
k2 = S.k2;
d = S.d;
I = S.I;
max_time = S.max_time;


% Initial conditions.
x0 = S.x0;
y0 = S.y0;
xdot0 = S.xdot0;
ydot0 = S.ydot0;

f0_bending = 1 / (2 * pi) * sqrt(k1 / (m1 + m2))
f0_torsion = 1 / (2 * pi) * sqrt(k2 / (I + m2 * d^2))

t = [0:0.0001:max_time];
x = nan(size(t));
u = nan(size(t));
udot = nan(size(t));
y = nan(size(t));
v = nan(size(t));
vdot = nan(size(t));
theta = nan(size(t));
omega = nan(size(t));
omegadot = nan(size(t));


x(1) = x0;
u(1) = xdot0;
udot(1) = 0;
y(1) = y0;
v(1) = ydot0;
vdot(1) = 0;
theta(1) = 0;
omega(1) = 0;
omegadot(1) = 0;


for i = 2 : length(t)
    dt = t(i) - t(i - 1);
    
    theta_now = theta(i - 1);
    omega_now = omega(i - 1);
    omegadot_now = omegadot(i - 1);
    x_now = x(i - 1);
    y_now = y(i - 1);
    
    % Differential equation for x
    udot(i) = 1 / (m1 + m2) * (cos(theta_now) * m2 * d * omegadot_now ...
        - sin(theta_now) * m2 * d * omega_now^2 ...
        - k1 * x_now);
    
    u(i) = u(i - 1) + udot(i) * dt;
    x(i) = x(i - 1) + u(i) * dt;
    
    % Differential equation for y
    vdot(i) = 1 / (m1 + m2) * (sin(theta_now) *  m2 * d * omegadot_now ...
        + cos(theta_now) * m2 * d * omega_now^2 ...
        - k1 * y_now);
    
    v(i) = v(i - 1) + vdot(i) * dt;
    y(i) = y(i - 1) + v(i) * dt;
    
    % Differential equation for theta
    omegadot(i) = 1 / (I + m2 * d^2) * (cos(theta_now) * m2 * d * udot(i) ...
        + sin(theta_now) * m2 * d * vdot(i) ...
        - damp_theta * omega_now ...
        - k2 * theta_now);
    
    omega(i) = omega(i - 1) + omegadot(i) * dt;
    theta(i) = theta(i - 1) + omega(i) * dt;
end

fig1 = figure('position', [100 100 1000 700]);
layout = tiledlayout(3, 4);
nexttile(1, [1, 4])
yyaxis left
h1 = plot(t, x, '--k', t, y, ':k');
ylabel('Position (m)')
ylim(max(abs(h1(1).Parent.YLim)).*[-1.01 1.01])
yyaxis right
h2 = plot(t, theta / pi * 180, '-b');
ylim(max(abs(h2.Parent.YLim)).*[-1 1])
ylabel('\theta (deg)');
xlabel('Time (s)');
legend('x', 'y', '\theta', 'box', 'off', 'orientation','horizontal', 'location', 'northoutside');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
box off

% Plot Lissajous figures.
nr_plots = 8;
fsize = 6;
start_time_indices = ceil([0 : nr_plots - 1] / nr_plots * length(t)) + 1;
end_time_indices = ceil([1 : nr_plots] / nr_plots * length(t));
for i = 1 : nr_plots
    nexttile()
    hold on
    plot(x, y, 'color', [1 1 1] * 0.7);
    start_time_index = start_time_indices(i);
    start_time = t(start_time_index);
    end_time_index = end_time_indices(i);
    end_time = t(end_time_index);
    plot(x(start_time_index:end_time_index), y(start_time_index:end_time_index), '-k', 'linewidth', 1.5);
    hold on
    %plot(x(start_time_index), y(start_time_index), 'xk', 'linewidth', 1)
    %xdelta = 0.08 * x(1);
    %text(x(start_time_index)+ xdelta, y(start_time_index), ['t = ' num2str(round(start_time, 1)) ' s'], 'color', 'k', 'fontsize', fsize)
    text(x(end_time_index)+ xdelta, y(end_time_index), ['t = ' num2str(end_time) ' s'], 'color', 'k', ...
        'backgroundcolor', 'white', 'fontsize', fsize)
    plot(x(end_time_index), y(end_time_index), 'ok', 'markerfacecolor', 'k');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal
    lim_pos = max(abs([min(x) max(x) min(y) max(y)])) * 1.1;
    xlim([-lim_pos lim_pos])
    ylim([-lim_pos lim_pos])
    box off
end

exportgraphics(layout, [S.label '.png'], 'resolution', 300)

if DO_ANIMATE
    fig3 = figure();
    hold on
    set(gcf,'color','w');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal
    
    % Thanks to: https://de.mathworks.com/matlabcentral/answers/385245-how-can-i-create-a-text-box-alongside-my-plot
    % set the width of the axis (the third value in Position) 
    % to be 60% of the Figure's width
    a = gca;
    a.Position(3) = 0.6;
    % put the textbox at 75% of the width and 
    % 10% of the height of the figure
    s = {['m_1 = ' num2str(m1) ' kg'], ...
        ['m_2 = ' num2str(m2) ' kg']...
        ['k_1 = ' num2str(k1) ' N/m']...
        ['k_2 = ' num2str(k2) ' Nm/rad']...
        ['d = ' num2str(d) ' m']...
        ['I_{zz} ' num2str(I) ' kg m^2']...
        };
    annotation('textbox', [0.75, 0.5, 0.1, 0.1], 'String', s)

    
    lim_val = sqrt(x0^2 + y0^2);
    factor = 1.1;
    axis([-lim_val, lim_val, -lim_val, lim_val] * factor)
    step_size = 50;
    gif_step_size = 1000;
    %h_line = animatedline;
    h_point = plot(x(i), y(i), 'ok', 'markerfacecolor', 'k');
    h_line = plot(x(i), y(i), '-k');
    for i = 1 : step_size : length(t)
        %addpoints(h_line, x(i), y(i));
        h_point.XData = x(i); %change x coordinate of the point
        h_point.YData = y(i); %change y coordinate of the point
        h_line.XData = [h_line.XData x(i)]; %add to line
        h_line.YData = [h_line.YData y(i)]; %add to line
        title([num2str(t(i)) ' s']);
        drawnow
        
        if DO_WRITE_GIF && mod(i, gif_step_size) == 1
            % Capture the plot as an image 
            frame = getframe(fig3); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,16); 
            % Write to the GIF File 
            if i == 1 
              imwrite(imind,cm,['trajectory_' S.label '.gif'],'gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,['trajectory_' S.label '.gif'],'gif','WriteMode','append'); 
            end 
        end
    end

end
