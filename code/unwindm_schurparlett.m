function UA=unwindm_schurparlett(A)
    %Creazione e opportuno riordinamento della forma di Schur.
    [U,T]=schur(A,'complex');
    [Uo,To]=ordschur(U,T,unwind(diag(T)));
    
    %Costruzione dell'array blocks contenente le dimensioni dei blocchi
    %nella forma di Schur riordinata.
    unwindings=unwind(diag(To));
    values=unique(unwindings,'stable');
    blocks=[];
    for i=1:length(values)
        blocks(i)=sum(unwindings==values(i));
    end
    
    %Costruzione dei blocchi diagonali.
    UTo=diag(unwindings);
    
    %Sottofunzione per l'indicizzazione dei blocchi.
    function r=rng(k)
        if k==1
            base_k=0;
        else
            base_k=sum(blocks(1:k-1));
        end
        r=(base_k+1):(base_k+blocks(k));
    end
    
    %Costruzione dei blocchi sovradiagonali.
    for j=2:length(blocks)
        for i=j-1:-1:1
            S1=To(rng(i),rng(i));
            S2=-To(rng(j),rng(j));
            S3=(values(i)-values(j))*To(rng(i),rng(j));
            for k=(i+1):(j-1)
                S3=S3+UTo(rng(i),rng(k))*To(rng(k),rng(j));
                S3=S3-To(rng(i),rng(k))*UTo(rng(k),rng(j));
            end
            UTo(rng(i),rng(j))=lyap(S1,S2,-S3);
        end
    end
    
    %Esecuzione della similitudine di Schur inversa.
    UA=Uo*UTo*Uo^(-1);
end
