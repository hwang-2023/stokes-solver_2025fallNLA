function [U, V, iter, time] = Inexact_uzawa(F, G, n, l, v1, v2, alpha, tau, MAXITER)

    tic;

    N = 2^n;
    r0 = sqrt(norm(F,'fro')^2 + norm(G,'fro')^2);

    U = zeros(N+1,N);
    V = zeros(N,N+1);
    P = zeros(N,N);

    rh = r0;
    
    iter = 0;

    epsilon = 1;

    while rh/r0 > 1e-8 && iter < MAXITER
        iter = iter + 1;

        [F1,G1] = left_multiply_B(P,N);

        [U,V,iter1] = PCG(F-F1,G-G1,N,2^l,v1,v2,epsilon);
        P = P + alpha*left_multiply_BT(U,V,N);
        
        epsilon = tau*norm(left_multiply_BT(U,V,N),'fro');

        [U1,V1,P1] = left_multiply_matrix(U,V,P,N);
        rh = sqrt(norm(U1-F,'fro')^2 + norm(V1-G,'fro')^2 + norm(P1,'fro')^2);

        disp(['Inexact Uzawa iteration',num2str(iter),': PCG iter=',num2str(iter1),', rh=',num2str(rh)]);
    end


    time = toc;
end
