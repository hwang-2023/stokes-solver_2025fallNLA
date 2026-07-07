MAXITER = 20;

time = zeros(1,4);
iter = zeros(1,4);
en = zeros(1,4);
logen = zeros(1,4);
N_list = [64, 128, 256, 512];


for k = 1:4
    n = k+5;
    N = 2^n;

    disp('Solve Stokes equation using Uzawa Iteration Method')
    disp(['N=',num2str(2^n),', alpha=',num2str(alpha),', MAXITER=',num2str(MAXITER)]);
    [F,G,U_true,V_true] = initialize(N);
    [U,V,iter(k),time(k)] = Uzawa(F,G,n,alpha,MAXITER);
    en(k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(k) = log(en(k));
    disp(['en=',num2str(en(k))]);
    disp(['iter=',num2str(iter(k)),', time=',num2str(time(k))]);

end


figure(1)
plot(N_list, time)
title('time')

figure(2)
plot(N_list, iter)
title('iter')

figure(3)
plot(N_list, en)
title('en')

figure(4)
plot(N_list, logen)
title('logen')

