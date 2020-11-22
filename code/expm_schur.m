function expA=expm_schur(A)
    [U,T]=schur(A,'complex');
    expA=U*expm(T)*U^(-1);
end

