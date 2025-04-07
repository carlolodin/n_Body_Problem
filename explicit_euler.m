%% Algoritmo Explicit Euler
function [newpos, newvel, newacc] = explicit_euler(oldpos, oldvel, dt, mass, G)
    newacc = compute_accelerations(oldpos, mass, G);
    newpos = oldpos + oldvel * dt;    
    newvel = oldvel + newacc * dt;
end

