function expA=expm_exact(A)
    if length(unique(eig(A)))~=length(eig(A))
        disp('attenzione, matrice non diagonalizzabile')
    end
    [Z,D]=eig(A);
    expD=diag(exp(diag(D)));
    expA=Z*expD*Z^(-1);
end
