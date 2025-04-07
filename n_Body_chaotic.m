% Simulazione del problema degli n corpi in 3D con velocity Verlet
clear; clc; close all;

%% Parametri
n = 4;                      % Numero di corpi
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

pos2 = pos;
vel2 = vel;

% Aggiunta di una piccola perturbazione
epsilon = 1e-4;
pos2(2:end,:) = pos2(2:end,:) + epsilon * randn(size(pos2(2:end,:)));


pos_hist2 = zeros(n,3,steps);
divergence = zeros(steps,1); % distanza media tra i due sistemi


%% Calcolo simulazione
acc = compute_accelerations(pos, mass, G);
acc2 = acc;
for t = 1:steps
    [pos, vel, acc] = verlet(pos, vel, acc, dt, mass, G); % USARE DT = 0.05
    % Salva posizione
    pos_hist(:,:,t) = pos;

    % Sistema 2
    [pos2, vel2, acc2] = verlet(pos2, vel2, acc2, dt, mass, G);
    pos_hist2(:,:,t) = pos2;

    % Calcola divergenza
    diff = pos - pos2;
    divergence(t) = mean(vecnorm(diff,2,2));  % media delle distanze tra corrispondenti corpi

    
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

figure('Name','Divergenza tra sistemi simili');
plot((1:steps)*dt, divergence, 'LineWidth', 1.5);
xlabel('Tempo');
ylabel('Distanza media tra i due sistemi');
title('Sensibilità alle condizioni iniziali (caos)');
grid on;


%% Plot animato in 3D con traiettorie che crescono nel tempo

figure('Name', 'Simulazione n-corpi');
colors = lines(n);
light_colors = colors + 0.5;             % schiarisce i colori
light_colors(light_colors > 1) = 1;      % mantieni entro [0,1]
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
trails2 = gobjects(n,1);
points2 = gobjects(n,1);
for i = 1:n
    % Sistema originale
    trails(i) = plot3(NaN,NaN,NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    points(i) = plot3(NaN,NaN,NaN, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 6);

    % Sistema perturbato (traiettoria chiara, punto trasparente)
    trails2(i) = plot3(NaN,NaN,NaN, '--', 'Color', light_colors(i,:), 'LineWidth', 1);
    points2(i) = plot3(NaN,NaN,NaN, 'o', 'MarkerFaceColor', light_colors(i,:), 'MarkerEdgeColor','none','MarkerSize', 4);
end


for t = 1:steps
    for i = 1:n
        % Aggiorna traiettoria progressiva originale
        set(trails(i), ...
            'XData', squeeze(pos_hist(i,1,1:t)), ...
            'YData', squeeze(pos_hist(i,2,1:t)), ...
            'ZData', squeeze(pos_hist(i,3,1:t)));

        % Punto originale
        set(points(i), ...
            'XData', pos_hist(i,1,t), ...
            'YData', pos_hist(i,2,t), ...
            'ZData', pos_hist(i,3,t));

        % Aggiorna traiettoria perturbata
        set(trails2(i), ...
            'XData', squeeze(pos_hist2(i,1,1:t)), ...
            'YData', squeeze(pos_hist2(i,2,1:t)), ...
            'ZData', squeeze(pos_hist2(i,3,1:t)));

        % Punto perturbato
        set(points2(i), ...
            'XData', pos_hist2(i,1,t), ...
            'YData', pos_hist2(i,2,t), ...
            'ZData', pos_hist2(i,3,t));
    end

    % ... (segue aggiornamento centro di massa come già nel tuo codice)

    drawnow;
end


