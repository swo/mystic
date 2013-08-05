s0 = test();

u0 = exp(s0(:,:,1));
u1 = exp(s0(:,:,3));

clf;
hold on;
surf(u0)
surf(u1)
hold off;