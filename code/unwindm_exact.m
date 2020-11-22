function UA=unwindm_exact(A)
    if length(unique(eig(A)))~=length(eig(A))
        disp('attenzione, matrice non diagonalizzabile')
    end
    [Z,D]=eig(A);
    UD=sym(diag(unwind(diag(D))));
    Z=sym(Z);
    UA=double(Z*UD*Z^(-1));
end