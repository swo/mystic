s = species_map();
[t, y] = run();

%to_show = {'O', 'C', 'N+', 'N-', 'S+', 'S-'};
to_show = {'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
idx = cellfun(@(x) s(x), to_show);

clf;

twod = 1;

hold all;
for i = idx
    plot(y(end, :, i))
end

xlabel('depth');
ylabel('time');
zlabel('concentration');

%legend('O', 'C', 'N+', 'N-', 'S+', 'S-');
legend('N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-');

hold off;