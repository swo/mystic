s = species_map();
sol = run();

clf;

l = 0;
hold on;
for i = [s('O(0)'), s('N(-III)'), s('C(IV)')]
    u = sol(:, :, i);
    
    if l
        u(u < 0) = 1e-7;
    end
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

if l
    set(gca, 'ZScale', 'log');
end

hold off;