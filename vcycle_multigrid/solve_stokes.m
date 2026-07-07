
n = 10;
l = 2;
v1 = 4;
v2 = 4;
MAXITER =20;

N = 2^n;

disp('Solve Stokes equation using V-cycle method')
disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);

[F, G, U_true, V_true] = initialize(2^n);

[U, V, iter, time] = V_cycle(F, G, v1, v2, n, l, MAXITER);


en = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);

disp(['en=',num2str(en)]);

disp(['iter=',num2str(iter),', time=',num2str(time)]);
