function expA=expm_mod_schur(A)
    [U,T]=schur(A,'complex');
    Tred=modm(T);
    if norm(Tred,'fro')>norm(T,'fro')
        disp('unwind(A) mal condizionata: è preferibile expm.')
    end
    expA=U*expm(Tred)*U^(-1);
end

