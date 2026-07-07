MAXITER = 20;

time = zeros(6,6);
iter = zeros(6,6);
en = zeros(6,6);
logen = zeros(6,6);
N_list = [64, 128, 256, 512, 1024, 2048];

for k = 1:6
    n = k+5;
    N = 2^n;

    v1=2;v2=2;l=1;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(1,k), time(1,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(1,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(1,k) = log(en(1,k));
    disp(['en=',num2str(en(1,k))]);
    disp(['iter=',num2str(iter(1,k)),', time=',num2str(time(1,k))]);

    v1=2;v2=2;l=2;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(2,k), time(2,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(2,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(2,k) = log(en(2,k));
    disp(['en=',num2str(en(2,k))]);
    disp(['iter=',num2str(iter(2,k)),', time=',num2str(time(2,k))]);

    v1=4;v2=4;l=1;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(3,k), time(3,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(3,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(3,k) = log(en(3,k));
    disp(['en=',num2str(en(3,k))]);
    disp(['iter=',num2str(iter(3,k)),', time=',num2str(time(3,k))]);

    v1=4;v2=4;l=2;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(4,k), time(4,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(4,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(4,k) = log(en(4,k));
    disp(['en=',num2str(en(4,k))]);
    disp(['iter=',num2str(iter(4,k)),', time=',num2str(time(4,k))]);

    v1=8;v2=8;l=1;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(5,k), time(5,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(5,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(5,k) = log(en(5,k));
    disp(['en=',num2str(en(5,k))]);
    disp(['iter=',num2str(iter(5,k)),', time=',num2str(time(5,k))]);

    v1=8;v2=8;l=2;

    disp('Solve Stokes equation using V-cycle method')
    disp(['N=',num2str(N),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', MAXITER=',num2str(MAXITER)]);
    [F, G, U_true, V_true] = initialize(2^n);
    [U, V, iter(6,k), time(6,k)] = V_cycle(F, G, v1, v2, n, l, MAXITER);
    en(6,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
    logen(6,k) = log(en(6,k));
    disp(['en=',num2str(en(6,k))]);
    disp(['iter=',num2str(iter(6,k)),', time=',num2str(time(6,k))]);

end


figure(1)
plot(N_list, time(1,:), N_list, time(2,:), N_list, time(3,:), N_list, time(4,:), N_list, time(5,:), N_list, time(6,:))
legend('v1=2,v2=2,l=1', 'v1=2,v2=2,l=2', 'v1=4,v2=4,l=1', 'v1=4,v2=4,l=2', 'v1=8,v2=8,l=1', 'v1=8,v2=8,l=2')
title('time')

figure(2)
plot(N_list, iter(1,:), N_list, iter(2,:), N_list, iter(3,:), N_list, iter(4,:), N_list, iter(5,:), N_list, iter(6,:))
legend('v1=2,v2=2,l=1', 'v1=2,v2=2,l=2', 'v1=4,v2=4,l=1', 'v1=4,v2=4,l=2', 'v1=8,v2=8,l=1', 'v1=8,v2=8,l=2')
title('iter')

figure(3)
plot(N_list, en(1,:), N_list, en(2,:), N_list, en(3,:), N_list, en(4,:), N_list, en(5,:), N_list, en(6,:))
legend('v1=2,v2=2,l=1', 'v1=2,v2=2,l=2', 'v1=4,v2=4,l=1', 'v1=4,v2=4,l=2', 'v1=8,v2=8,l=1', 'v1=8,v2=8,l=2')
title('en')

figure(4)
plot(N_list, logen(1,:), N_list, logen(2,:), N_list, logen(3,:), N_list, logen(4,:), N_list, logen(5,:), N_list, logen(6,:))
legend('v1=2,v2=2,l=1', 'v1=2,v2=2,l=2', 'v1=4,v2=4,l=1', 'v1=4,v2=4,l=2', 'v1=8,v2=8,l=1', 'v1=8,v2=8,l=2')
title('logen')

disp('time');
disp('v1=2,v2=2,l=1');
disp(time(1,:));
disp('v1=2,v2=2,l=2');
disp(time(2,:));
disp('v1=4,v2=4,l=1');
disp(time(3,:));
disp('v1=4,v2=4,l=2');
disp(time(4,:));
disp('v1=8,v2=8,l=1');
disp(time(5,:));
disp('v1=8,v2=8,l=2');
disp(time(6,:));

disp('iter');
disp('v1=2,v2=2,l=1');
disp(iter(1,:));
disp('v1=2,v2=2,l=2');
disp(iter(2,:));
disp('v1=4,v2=4,l=1');
disp(iter(3,:));
disp('v1=4,v2=4,l=2');
disp(iter(4,:));
disp('v1=8,v2=8,l=1');
disp(iter(5,:));
disp('v1=8,v2=8,l=2');
disp(iter(6,:));

disp('en');
disp('v1=2,v2=2,l=1');
disp(en(1,:));
disp('v1=2,v2=2,l=2');
disp(en(2,:));
disp('v1=4,v2=4,l=1');
disp(en(3,:));
disp('v1=4,v2=4,l=2');
disp(en(4,:));
disp('v1=8,v2=8,l=1');
disp(en(5,:));
disp('v1=8,v2=8,l=2');
disp(en(6,:));

disp('logen');
disp('v1=2,v2=2,l=1');
disp(logen(1,:));
disp('v1=2,v2=2,l=2');
disp(logen(2,:));
disp('v1=4,v2=4,l=1');
disp(logen(3,:));
disp('v1=4,v2=4,l=2');
disp(logen(4,:));
disp('v1=8,v2=8,l=1');
disp(logen(5,:));
disp('v1=8,v2=8,l=2');
disp(logen(6,:));

