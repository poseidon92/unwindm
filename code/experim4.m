X=diag(100i*[1:10]);
U=gallery('randsvd',10,500);
B=U*X*U^(-1);

scalings=[];
scalings_mod=[];
err=[];
err_mod=[];
found=0;
t=0:1:100;
t(1)=1;
for i=1:length(t)
    [UB,TB]=schur(B,'complex');
    scalings(i)=scalingfactor(t(i)*TB);
    scalings_mod(i)=scalingfactor(t(i)*(modm(TB)));
    exact=expm_exact(t(i)*B);
    expB=expm_schur(t(i)*B);
    expBmod=expm_mod_schur(t(i)*B);
    err(i)=norm(expB-exact,2)/norm(exact,2);
    err_mod(i)=norm(expBmod-exact,2)/norm(exact,2);
    if found==0
        if scalings_mod(i)<scalings(i)-4
            found=1;
            gainpos=i;
        end
    end
end

figure(1);
title('Scaling')
hold on
xline(t(gainpos), 'LineStyle','--');
plot(t,scalings, 'Color', '#3743f7');
plot(t,scalings_mod, 'Color', '#809fb7');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('fattore di scaling $s$', 'Interpreter', 'latex');
lg=legend('soglia di convenienza computazionale','algoritmo (a)','algoritmo (b)');
set(lg,'Interpreter','latex','Location','southeast');
exportgraphics(gca,'experim4_1.png','Resolution',300);

figure(2);
semilogy(t,err,'Color', '#3743f7','Marker','o','MarkerSize',12,'LineStyle','none');
title('Errori relativi');
hold on
semilogy(t,err_mod,'Color','#809fb7','Marker','x','MarkerSize',12,'LineStyle','none');
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\epsilon_{rel}$','Interpreter','latex')
lg=legend('$\epsilon_{rel}$ in alg. (a)','$\epsilon_{rel}$ in alg. (b)');
set(lg,'Interpreter','latex','Location','southeast');
exportgraphics(gca,'experim4_2.png','Resolution',300);

