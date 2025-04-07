%% Algoritmo Runge Kutta 4
function [newpos, newvel, newacc] = runge_kutta_four(oldpos, oldvel, dt, mass, G)
    n = length(mass);
    y = [oldpos(:); oldvel(:)];
    k1 = derivatives(y, mass, G);
    k2 = derivatives(y + 0.5*dt*k1, mass, G);
    k3 = derivatives(y + 0.5*dt*k2, mass, G);
    k4 = derivatives(y + dt*k3, mass, G);
    y = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    newpos = reshape(y(1:n*3), [n,3]);
    newvel = reshape(y(n*3+1:end), [n,3]);
    newacc = compute_accelerations(newpos, mass, G);
end


