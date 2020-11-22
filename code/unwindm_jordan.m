function [UA,A_sym,Z,J]=unwindm_jordan(A)
    A_sym=sym(A); 
    [Z,J]=jordan(A_sym);
    UA=Z*diag(unwind(diag(J)))*Z^(-1);
end

