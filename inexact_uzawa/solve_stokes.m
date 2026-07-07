n = 10;
l = 1;
N = 2^n;
L = 2^l; 
v1 = 4;
v2 = 4;
alpha = 1;
tau = 1e-3;
MAXITER = 30;


disp('Solve Stokes equation using Inexact Uzawa Iteration Method')
disp(['N=',num2str(2^n),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', alpha=',num2str(alpha),', tau=',num2str(tau),', MAXITER=',num2str(MAXITER)]);

[F,G,U_true,V_true] = initialize(N);

[U,V,iter,time] = Inexact_uzawa(F,G,n,l,v1,v2,alpha,tau,MAXITER);


en = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);

disp(['en=',num2str(en)]);

disp(['iter=',num2str(iter),', time=',num2str(time)]);

