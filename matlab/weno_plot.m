
figure, hold on;
plot(weno(:,1), 'ko-')
plot(weno(:,2), 'b*-')
plot(weno(:,3), 'r')
legend('function', 'analytic', 'weno', 'location', 'southeast')