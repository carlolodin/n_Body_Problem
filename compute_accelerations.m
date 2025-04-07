%% Funzione per calcolare accelerazioni
function acc = compute_accelerations(pos, mass, G, t)
    n = size(pos,1);
    acc = zeros(n,3);
    for i = 1:n
        for j = 1:n
            if i ~= j
                r = pos(j,:) - pos(i,:);
                dist = norm(r) + 1e-5; % Evita divisione per zero
                acc(i,:) = acc(i,:) + G * mass(j) * r / dist^3;
            end
        end
    end
end

