A=i*[-1/2 -7/2 17/2;-6 2 9;1/2 -1/2 13/2];
[UA,A_sym,Z,J]=unwindm_jordan(A)

B=[0 5 0;-5 0 5;0 -5 0];
[UB,B_sym,Z,J]=unwindm_jordan(B)