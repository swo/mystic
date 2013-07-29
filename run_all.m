sol = run();

clf;

hold on;
for i = [10]
    surf(sol(:, :, i))
end
xlabel('depth');
ylabel('time');
zlabel('concentration');
%set(gca, 'ZScale', 'log');
hold off;