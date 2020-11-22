scalings=[];
scalings_mod=[];
err=[];
err_mod=[];
fails=[];
found=0;
k=1:1:20;
for i=1:length(k)
    X=convdiffmatrix(50,k(i));
    [U,T]=schur(X,'complex');
    scalings(i)=scalingfactor(T);
    scalings_mod(i)=scalingfactor(modm(T));
    if scalings_mod(i) > scalings(i)
        fails=cat(1,fails,i);
    end
    exact=expm_exact(X);
    expX=expm_schur(X);
    expXmod=expm_mod_schur(X);
    err(i)=norm(expX-exact,2)/norm(exact,2);
    err_mod(i)=norm(expXmod-exact,2)/norm(exact,2);
    if found==0
        if scalings_mod(i)<scalings(i)-4
            found=1;
            gainpos=i;
        end
    end
end

figure(1);
title('Convection-diffusion: scaling')
hold on
xline(k(gainpos), 'LineStyle','--');
plot(k,scalings, 'Color', '#3743f7');
plot(k,scalings_mod, 'Color', '#809fb7');
plot(fails,scalings_mod(fails),'LineStyle','none','Marker','x','MarkerSize',12,'Color','#cc123b')
xlabel('$k$ tale che $c=1.6^k$ e $d=0.2(0.5)^k$', 'Interpreter', 'latex');
ylabel('fattore di scaling $s$', 'Interpreter', 'latex');
lg=legend('soglia di convenienza computazionale','algoritmo (a)','algoritmo (b)','breakdown di (b)');
set(lg,'Interpreter','latex','Location','southeast');
exportgraphics(gca,'experim5_1.png','Resolution',600);

figure(2);
st=max(fails)+1;
semilogy(k,err,'Color', '#3743f7','Marker','o','MarkerSize',12,'LineStyle','none');
title('Convection-diffusion: errori relativi');
hold on
semilogy(k(st:end),err_mod(st:end),'Color','#809fb7','Marker','x','MarkerSize',12,'LineStyle','none');
xlabel('$k$ tale che $c=1.6^k$ e $d=0.2(0.5)^k$', 'Interpreter', 'latex')
ylabel('$\epsilon_{rel}$','Interpreter','latex')
lg=legend('$\epsilon_{rel}$ in alg. (a)','$\epsilon_{rel}$ in alg. (b)');
set(lg,'Interpreter','latex','Location','southeast');
exportgraphics(gca,'experim5_2.png','Resolution',300);

