load('tols1090.mat');
hold on
title('Spettro della matrice Tolosa1090');
plot(eig(tols1090),'Color','#3743f7','LineStyle','none','Marker','o','MarkerSize',12);
xlabel('$\Re(z)$','Interpreter','latex');
ylabel('$\Im(z)$','Interpreter','latex');
exportgraphics(gca,'tols1090spectrum.png','Resolution',600); 