%% Funzione derivatives
function dy = derivatives(y, mass, G)
    n = length(mass);
    pos = reshape(y(1:n*3), [n, 3]);
    vel = reshape(y(n*3+1:end), [n, 3]);
    acc = compute_accelerations(pos, mass, G);
    dy = [vel(:); acc(:)];
end