% Simulazione del problema degli n corpi in 3D con velocity Verlet
clear; clc; close all;

%% Parametri
n = 10;                      % Numero di corpi
G = 1;                      % Costante gravitazionale (normalizzata)
T = 20;                     % Tempo totale della simulazione
dt = 0.01;                  % Passo temporale
steps = floor(T/dt);        % Numero di passi

cm_pos = zeros(steps, 3);
cm_vel = zeros(steps, 3);
L = zeros(steps, 1);        % Modulo del momento angolare totale
v_cm_mag = zeros(steps, 1); % Modulo velocità centro di massa
pos_hist = zeros(n,3,steps);
energy = zeros(steps,1);


mass = [1; ones(n-1, 1)];

% Corpo centrale al centro
pos = zeros(n,3);
vel = zeros(n,3);

% Distribuisci gli altri corpi su una sfera attorno al centro
radii = linspace(2, 5, n-1)';
phi = rand(n-1,1) * 2*pi;         % angolo longitudinale
theta = acos(2*rand(n-1,1) - 1);  % angolo latitudinale (uniforme sulla sfera)

% Conversione coordinate sferiche → cartesiane
pos(2:end,1) = radii .* sin(theta) .* cos(phi);
pos(2:end,2) = radii .* sin(theta) .* sin(phi);
pos(2:end,3) = radii .* cos(theta);

% Assegna velocità tangenziali per creare orbite 3D
for i = 2:n
    r_vec = pos(i,:) - pos(1,:);
    r = norm(r_vec);
    v_mag = sqrt(G * mass(1) / r);

    % Calcola un vettore perpendicolare a r_vec
    if abs(dot(r_vec, [0 0 1])) < 0.9
        temp = [0 0 1];
    else
        temp = [0 1 0];
    end
    tangent = cross(r_vec, temp);
    tangent = tangent / norm(tangent);

    vel(i,:) = v_mag * tangent;
end

% Il corpo centrale è fermo
vel(1,:) = [0 0 0];



%% Calcolo simulazione
acc = compute_accelerations(pos, mass, G);
for t = 1:steps
    % Choose a numerical method
    [pos, vel, acc] = runge_kutta_four(pos, vel, dt, mass, G);
    %[pos, vel, acc] = verlet(pos, vel, acc, dt, mass, G); % USARE DT = 0.05
    %[pos, vel, acc] = symplectic_euler(pos, vel, dt, mass, G);
    %[pos, vel, acc] = explicit_euler(pos, vel, dt, mass, G);

    % Salva posizione
    pos_hist(:,:,t) = pos;
    
    % Calcola energia totale
    KE = 0.5 * sum(mass .* sum(vel.^2, 2));
    PE = 0;
    for i = 1:n
        for j = i+1:n
            r = norm(pos(i,:) - pos(j,:));
            PE = PE - G * mass(i) * mass(j) / r;
        end
    end
    energy(t) = KE + PE;

    % Calcolo del centro di massa
    cm_pos(t,:) = sum(mass .* pos) / sum(mass);
    cm_vel(t,:) = sum(mass .* vel) / sum(mass);

    % Modulo velocità CM
    v_cm_mag(t) = norm(cm_vel(t,:));

    % Momento angolare totale
    L_vec = zeros(1,3);
    for i = 1:n
        L_vec = L_vec + mass(i) * cross(pos(i,:), vel(i,:));
    end
    L(t) = norm(L_vec);

end

%% Plot energia, momento angolare e velocità del CM
figure('Name', 'Conservazioni');

t_vec = (1:steps)*dt;

subplot(3,1,1);
plot(t_vec, energy, 'b-', 'LineWidth', 1.5);
xlabel('Tempo');
ylabel('Energia Totale');
title('Energia Totale del Sistema');
grid on;

subplot(3,1,2);
plot(t_vec, L, 'r-', 'LineWidth', 1.5);
xlabel('Tempo');
ylabel('|L|');
title('Modulo del Momento Angolare Totale');
grid on;

subplot(3,1,3);
plot(t_vec, v_cm_mag, 'k-', 'LineWidth', 1.5);
xlabel('Tempo');
ylabel('|V_{CM}|');
title('Modulo della Velocità del Centro di Massa');
grid on;

%% Plot animato in 3D con traiettorie che crescono nel tempo

figure('Name', 'Simulazione n-corpi');
colors = lines(n);
hold on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
view(3);

% Calcola i limiti globali per tutta la simulazione
all_pos = reshape(pos_hist, n, 3, []);
min_vals = min(min(all_pos,[],3),[],1);
max_vals = max(max(all_pos,[],3),[],1);
margin = 1;
lims = [min_vals - margin; max_vals + margin];
xlim([lims(1,1), lims(2,1)]);
ylim([lims(1,2), lims(2,2)]);
zlim([lims(1,3), lims(2,3)]);
axis manual;

% Tracce (che si allungano) e punti
trails = gobjects(n,1);
points = gobjects(n,1);
for i = 1:n
    trails(i) = plot3(NaN,NaN,NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    points(i) = plot3(NaN,NaN,NaN, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 6);
end

% Traccia del centro di massa
cm_trail = plot3(NaN, NaN, NaN, '--', 'Color', [1 0.6 0], 'LineWidth', 1.5);

% Punto e vettore velocità del centro di massa
cm_point = plot3(NaN, NaN, NaN, 'p', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.6 1], 'MarkerEdgeColor', 'k');

scale = 2;
cm_quiver = quiver3(NaN, NaN, NaN, NaN, NaN, NaN, 0, 'Color', [1 0.5 0], 'LineWidth', 2); % arancione

% Animazione
for t = 1:steps
    for i = 1:n
        % Aggiorna traiettoria progressiva
        set(trails(i), ...
            'XData', squeeze(pos_hist(i,1,1:t)), ...
            'YData', squeeze(pos_hist(i,2,1:t)), ...
            'ZData', squeeze(pos_hist(i,3,1:t)));

        % Aggiorna punto
        set(points(i), ...
            'XData', pos_hist(i,1,t), ...
            'YData', pos_hist(i,2,t), ...
            'ZData', pos_hist(i,3,t));
    end

    % Aggiorna CM
    set(cm_point, ...
        'XData', cm_pos(t,1), ...
        'YData', cm_pos(t,2), ...
        'ZData', cm_pos(t,3));

    % Traccia CM progressiva
    set(cm_trail, ...
        'XData', cm_pos(1:t,1), ...
        'YData', cm_pos(1:t,2), ...
        'ZData', cm_pos(1:t,3));

    % Vettore velocità CM
    set(cm_quiver, ...
        'XData', cm_pos(t,1), ...
        'YData', cm_pos(t,2), ...
        'ZData', cm_pos(t,3), ...
        'UData', cm_vel(t,1)*scale, ...
        'VData', cm_vel(t,2)*scale, ...
        'WData', cm_vel(t,3)*scale);

    drawnow;
end