% Simulazione del problema degli n corpi in 3D con velocity Verlet
clear; clc; close all;

%% Parametri
n = 5;                      % Numero di corpi
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


mass = [5; ones(n-1, 1)];

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

% Salva le condizioni iniziali per ogni metodo
pos0 = pos;
vel0 = vel;
acc0 = compute_accelerations(pos0, mass, G);


%% Calcolo simulazione
methods = {'EulerEsplicito', 'EulerSimplettico', 'Verlet', 'RungeKutta4'};
results = struct();

for m = 1:length(methods)
    method = methods{m};

    pos = pos0;
    vel = vel0;
    acc = acc0;

    % Prealloca
    pos_hist = zeros(n,3,steps);
    energy = zeros(steps,1);
    L = zeros(steps,1);
    v_cm_mag = zeros(steps,1);

    for t = 1:steps
        switch method
            case 'EulerEsplicito'
                [pos, vel, acc] = explicit_euler(pos, vel, dt, mass, G);
            case 'EulerSimplettico'
                [pos, vel, acc] = symplectic_euler(pos, vel, dt, mass, G);
            case 'Verlet'
                [pos, vel, acc] = verlet(pos, vel, acc, dt, mass, G);
            case 'RungeKutta4'
                [pos, vel, acc] = runge_kutta_four(pos, vel, dt, mass, G);
        end

        % Salva posizione
        pos_hist(:,:,t) = pos;

        % Energia, CM, Momento angolare (come prima)
        KE = 0.5 * sum(mass .* sum(vel.^2, 2));
        PE = 0;
        for i = 1:n
            for j = i+1:n
                r = norm(pos(i,:) - pos(j,:));
                PE = PE - G * mass(i) * mass(j) / max(r, 1e-4);
            end
        end
        energy(t) = KE + PE;

        cm_vel = sum(mass .* vel) / sum(mass);
        v_cm_mag(t) = norm(cm_vel);

        L_vec = zeros(1,3);
        for i = 1:n
            L_vec = L_vec + mass(i) * cross(pos(i,:), vel(i,:));
        end
        L(t) = norm(L_vec);
    end

    % Salva i risultati per questo metodo
    results.(method).energy = energy;
    results.(method).L = L;
    results.(method).v_cm_mag = v_cm_mag;
    results.(method).pos_hist = pos_hist;
end


%% Plot energia, momento angolare e velocità del CM
figure;
hold on;
for m = 1:length(methods)
    plot((1:steps)*dt, results.(methods{m}).energy, 'DisplayName', methods{m});
end
title('Energia Totale');
xlabel('Tempo'); ylabel('Energia');
legend; grid on;
%% Plot animato in 3D con traiettorie che crescono nel tempo

figure;
for m = 1:length(methods)
    subplot(2,2,m);
    hold on;
    for i = 1:n
        plot3(squeeze(results.(methods{m}).pos_hist(i,1,:)), ...
              squeeze(results.(methods{m}).pos_hist(i,2,:)), ...
              squeeze(results.(methods{m}).pos_hist(i,3,:)), ...
              'DisplayName', sprintf('Body %d', i));
    end
    title(['Traiettorie - ', methods{m}]);
    axis equal;
    grid on;
    view(3);
end