function [M,kfin]=badcondmatrix(N,k,pert)
    eigvals=(2*randi([0,5],N,1)+1)*i*pi+pert*(randn(N,1)+i*randn(N,1));
    Z=gallery('randsvd',N,k);
    Zinv=Z^(-1);
    D=diag(eigvals);
    [X,M]=schur(Z*D*Zinv,'complex');
    kfin=norm(Z,2)*norm(Zinv,2);
end

