n = 9;

N = 2^n;
alpha = 1;
MAXITER = 20;

disp('Solve Stokes equation using Uzawa Iteration Method')
disp(['N=',num2str(2^n),', alpha=',num2str(alpha),', MAXITER=',num2str(MAXITER)]);


[F,G,U_true,V_true] = initialize(N);

[U,V,iter,time] = Uzawa(F,G,n,alpha,MAXITER);


en = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);

disp(['en=',num2str(en)]);

disp(['iter=',num2str(iter),', time=',num2str(time)]);


