A=[1 -500;
   500 1];
B=[0 30 1 1 1 1;
   -100 0 1 1 1 1;
   0 0 0 -6 1 1;
   0 0 500 0 1 1;
   0 0 0 0 0 200;
   0 0 0 0 -15 0];
load('tols1090.mat');

[UA,TA]=schur(A,'complex');
[UB,TB]=schur(B,'complex');
[Utols,Ttols]=schur(tols1090,'complex');


resA=[];
resB=[];
restols1090=[];
t=[1,50,100];
for i=1:length(t)
    resA(i,1)=scalingfactor(i*TA);
    resA(i,2)=scalingfactor(i*(modm(TA)));
    resB(i,1)=scalingfactor(i*TB);
    resB(i,2)=scalingfactor(i*(modm(TB)));
    restols1090(i,1)=scalingfactor(i*Ttols);
    restols1090(i,2)=scalingfactor(i*(modm(Ttols)));
end

results.resA=resA;
results.resB=resB;
results.restols1090=restols1090;
save('experim2.mat','results');
    