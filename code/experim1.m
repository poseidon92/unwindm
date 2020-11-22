klist=[2,10,50,100];
rep=25;
trials=[];
for l=1:length(klist)
    trials=cat(1,trials,klist(l)*ones(rep,1));
end

erel_sp=[];
erel_le=[];
cond_lb=[];
for l=1:length(trials)
    lb=0;
    while lb<1e-12
        [A,k]=badcondmatrix(10,trials(l),1e-6);
        UA=unwindm_exact(A);
        UA_schurparlett=unwindm_schurparlett(A);
        UA_logexp=unwindm_logexp(A);
        erel_sp(l)=norm(UA_schurparlett-UA,2)/norm(UA,2);
        erel_le(l)=norm(UA_logexp-UA,2)/norm(UA,2);
        cond_lb(l)=cond_unwindm_lb(A,k);
        lb=cond_lb(l)*eps;
    end
end
max(max(erel_sp./cond_lb),max(erel_le./cond_lb)) %1.053888555483874e-13
figure(1);
semilogy(erel_sp,'Color','red','Marker','x','LineStyle','none','MarkerSize',12);
hold on
ylim([1e-16 1e-6])
semilogy(erel_le,'Color','blue','Marker','o','LineStyle','none','MarkerSize',12);
semilogy(cond_lb*eps,'Color','black','LineStyle','-.');
lg=legend('$\epsilon_{rel}$ Schur-Parlett','$\epsilon_{rel}$ log-exp','$cond^{-}_{\mathcal{U}} \cdot u$');
set(lg,'Interpreter','latex','Location','northwest');
exportgraphics(gca,'experim1.png','Resolution',300);  