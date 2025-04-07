%% Algoritmo verlet
function [newpos, newvel, newacc] = verlet(oldpos, oldvel, oldacc, dt, mass, G)
    % USARE DT = 0.05
    
    % Aggiorna posizioni
    newpos = oldpos + oldvel*dt + 0.5*oldacc*dt^2;
    
    % Calcola nuove accelerazioni
    newacc = compute_accelerations(newpos, mass, G);
    
    % Aggiorna velocità
    newvel = oldvel + 0.5*(oldacc + newacc)*dt;
    

end

