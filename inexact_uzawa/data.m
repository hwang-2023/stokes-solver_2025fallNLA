MAXITER = 30;
alpha=1;
tau=1e-3;


time = zeros(6,6);
iter = zeros(6,6);
en = zeros(6,6);
logen = zeros(6,6);
N_list = [64, 128, 256, 512, 1024, 2048];
l_list = [1,2];
v_list = [2,4];
L_list = [2,4];
alpha_list = [0.95,1,1.05];
tau_list = [1e-3,1e-4];

s = cell(1,24);
cnt = 0;
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                cnt = cnt + 1;
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                l = l_list(i3);                
                tau = tau_list(i4);
                for k = 1:6
                    n = k+5;
                    N = 2^n;
                    disp('Solve Stokes equation using Inexact Uzawa Iteration Method')
                    s{cnt} = sprintf('L=%d,v1=%d,v2=%d,alpha=%.2f,tau=%.4f',l,v1,v2,alpha,tau);
                    disp(['N=',num2str(2^n),', L=',num2str(2^l),', v1=',num2str(v1),', v2=',num2str(v2),', alpha=',num2str(alpha),', tau=',num2str(tau),', MAXITER=',num2str(MAXITER)]);
                    [F,G,U_true,V_true] = initialize(N);
                    [U,V,iter(j,k),time(j,k)] = Inexact_uzawa(F,G,n,l,v1,v2,alpha,tau,MAXITER);
                    en(j,k) = 1/N*sqrt(norm(U-U_true,'fro')^2+norm(V-V_true,'fro')^2);
                    logen(j,k) = log(en(j,k));
                    disp(['en=',num2str(en(j,k))]);
                    disp(['iter=',num2str(iter(j,k)),', time=',num2str(time(j,k))]);
                end
            end
        end
    end
end


figure(1);
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                plot(N_list, time(j,:))               
                hold on;
            end
        end
    end
end
title('time');
legend(s,'FontSize',7);

figure(2);
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                plot(N_list, iter(j,:))               
                hold on;
            end
        end
    end
end
title('iter');
legend(s,'FontSize',7);

figure(3);
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                plot(N_list, en(j,:))               
                hold on;
            end
        end
    end
end
title('en');
legend(s,'FontSize',7);

figure(4);
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                plot(N_list, logen(j,:))               
                hold on;
            end
        end
    end
end
title('logen');
legend(s,'FontSize',7);


cnt = 0;
disp('time');
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                cnt = cnt + 1;
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                disp(s{cnt});
                disp(time(j,:));
            end
        end
    end
end

cnt = 0;
disp('iter');
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                cnt = cnt + 1;
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                disp(s{cnt});
                disp(iter(j,:));
            end
        end
    end
end


cnt = 0;
disp('en');
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                cnt = cnt + 1;
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                disp(s{cnt});
                disp(en(j,:));
            end
        end
    end
end

cnt = 0;
disp('logen');
for i1 = 1:3
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                cnt = cnt + 1;
                j = 8*(i1-1) + 4*(i2-1) + 2*(i3-1) + (i4-1) + 1;
                alpha = alpha_list(i1);
                v1 = v_list(i2);
                v2 = v_list(i2);
                L = L_list(i3);                
                tau = tau_list(i4);
                disp(s{cnt});
                disp(logen(j,:));
            end
        end
    end
end