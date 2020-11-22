function M=convdiffmatrix(N,k)
    c=1.6^k;
    d=0.2*(0.5)^k;
    M=-2*d*eye(N) + diag((d+c)*ones(N-1,1),-1) + diag((d-c)*ones(N-1,1),1);
end
