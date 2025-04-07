%% Algoritmo Symplectic Euler
function [newpos, newvel, newacc] = symplectic_euler(oldpos, oldvel, dt, mass, G)
    newacc = compute_accelerations(oldpos, mass, G);
    newvel = oldvel + newacc * dt;
    newpos = oldpos + newvel * dt;
end

