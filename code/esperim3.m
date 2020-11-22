load('tols1090.mat')
scalings=[];
scalings_mod=[];
found=0;
t=0:1:30;
t(1)=1;
for i=1:length(t)
    [U,T]=schur(tols1090,'complex');
    scalings(i)=scalingfactor(t(i)*T);
    scalings_mod(i)=scalingfactor(t(i)*(modm(T)));
    if found==0
        if scalings_mod(i)<scalings(i)-4
            found=1;
            gainpos=i;
        end
    end
end

figure(1);
title('Matrice Tolosa1090: scaling')
hold on
xline(t(gainpos), 'LineStyle','--');
plot(t,scalings, 'Color', '#3743f7');
plot(t,scalings_mod, 'Color', '#809fb7');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('fattore di scaling $s$', 'Interpreter', 'latex');
lg=legend('soglia di convenienza computazionale','algoritmo (a)','algoritmo (b)');
set(lg,'Interpreter','latex','Location','southeast');
exportgraphics(gca,'experim3.png','Resolution',300);