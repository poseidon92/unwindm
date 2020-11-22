function s=scalingfactor(M)
    theta13=5.4;
    s=max(ceil(log2(norm(M,1)/theta13)),0);
end

