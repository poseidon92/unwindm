function cond=cond_unwindm_lb(A,k)
    %Calcolo il massimo modulo delle divided differences.
    mxabs=0;
    eigvls=eig(A);
    for i=1:length(eigvls)
        for j=i+1:length(eigvls)
            if eigvls(i)~=eigvls(j)
                mxabs=max(mxabs,abs((unwind(eigvls(i))-unwind(eigvls(j)))/(eigvls(i)-eigvls(j))));
            end
        end
    end
    
    cond=pi/k*mxabs;
end

