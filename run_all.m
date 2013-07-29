sol = run();

clf;

hold on;
for i = [2 3 6 7]
    surf(sol(:, :, i))
end
xlabel('depth');
ylabel('time');
zlabel('concentration');
%set(gca, 'ZScale', 'log');
hold off;