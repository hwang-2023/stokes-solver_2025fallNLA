function [U_half,V_half,P]=DGS(U,V,P,F,G,Z,N)

%实现DGS迭代
%U,V,P:初值
%F,G,Z:求解的方程的RHS

    h = 1/N;
   
    
    [U1,V1,~] = left_multiply_matrix(U,V,P,N); 

    %更新U

    %计算残差r
    RF=F-U1;


    %计算E=U'-U
    EU = zeros(N+1,N);

    EU(1,1) = RF(1,1);
    for i=2:N
        EU(i,1) = 1/3*EU(i-1,1) + 1/3*h^2*RF(i,1);
    end
    EU(N+1,1) = RF(N+1,1);

    for j = 2:N-1
        EU(1,j) = RF(1,j);
        for i = 2:N
            EU(i,j) = 1/4*EU(i,j-1) + 1/4*EU(i-1,j) + 1/4*h^2*RF(i,j);
        end
        EU(N+1,j) = RF(N+1,j);
    end
    
    EU(1,N) = RF(1,N);
    for i = 2:N
        EU(i,N) = 1/3*EU(i-1,N) + 1/3*EU(i,N-1) + 1/3*h^2*RF(i,N);
    end
    EU(N+1,1) = RF(N+1,1);

    %得U_half
    U_half = U + EU;
    

    %更新V
    RG=G-V1;
    
    
    EV = zeros(N,N+1);


    EV(1,1) = RG(1,1);
    for j=2:N
        EV(1,j) = 1/3*EV(1,j-1) + 1/3*h^2*RG(1,j);
    end
    EV(1,N+1) = RG(1,N+1);

    
    
    for i=2:N-1
        EV(i,1) = RG(i,1);
        for j=2:N
            EV(i,j) = 1/4*EV(i,j-1) + 1/4*EV(i-1,j) + 1/4*h^2*RG(i,j);
        end
        EV(i,N+1) = RG(i,N+1);
    end

    EV(N,1) = RG(N,1);
    for j=2:N
        EV(N,j) = 1/3*EV(N,j-1) + 1/3*h^2*RG(N,j);
    end
    EV(N,N+1) = RG(N,N+1);

    
    
    V_half = V + EV;


    %更新P

    %内部
    for i = 2:N-1
        for j = 2:N-1
            r = -1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j)) - Z(i,j);
            delta = r*h/4;
            U_half(i,j) = U_half(i,j) - delta;
            U_half(i+1,j) = U_half(i+1,j) + delta;
            V_half(i,j) = V_half(i,j) - delta;
            V_half(i,j+1) = V_half(i,j+1) + delta;
            P(i,j) = P(i,j) + r;
            P(i+1,j) = P(i+1,j) - 1/4*r;
            P(i-1,j) = P(i-1,j) - 1/4*r;
            P(i,j+1) = P(i,j+1) - 1/4*r;
            P(i,j-1) = P(i,j-1) - 1/4*r;
        end 
    end

    %顶边与底边
    
    for i = 2:N-1
        r = -1/h*(U_half(i+1,1)-U_half(i,1)+V_half(i,2)-V_half(i,1)) - Z(i,1);
        delta = r*h/3;
        U_half(i,1) = U_half(i,1) - delta;
        U_half(i+1,1) = U_half(i+1,1) + delta;
        V_half(i,2) = V_half(i,2) + delta;
        P(i,1) = P(i,1) + r;
        P(i+1,1) = P(i+1,1) - 1/3*r;
        P(i-1,1) = P(i-1,1) - 1/3*r;
        P(i,2) = P(i,2) - 1/3*r;

        r = -1/h*(U_half(i+1,N)-U_half(i,N)+V_half(i,N+1)-V_half(i,N)) - Z(i,N);
        delta = r*h/3;
        U_half(i,N) = U_half(i,N) - delta;
        U_half(i+1,N) = U_half(i+1,N) + delta;
        V_half(i,N) = V_half(i,N) - delta;
        P(i,N) = P(i,N) + r;
        P(i+1,N) = P(i+1,N) - 1/3*r;
        P(i-1,N) = P(i-1,N) - 1/3*r;
        P(i,N-1) = P(i,N-1) - 1/3*r;
        
        
    end

    %左边与右边

    for j=2:N-1    
        r = -1/h*(U_half(2,j)-U_half(1,j)+V_half(1,j+1)-V_half(1,j)) - Z(1,j);
        delta = r*h/3;
        V_half(1,j) = V_half(1,j) - delta;
        V_half(1,j+1) = V_half(1,j+1) + delta;
        U_half(2,j) = U_half(1+1,j) + delta;
        P(1,j) = P(1,j) + r;
        P(1,j+1) = P(1,j+1) - 1/3*r;
        P(2,j) = P(2,j) - 1/3*r;
        P(1,j-1) = P(1,j-1) - 1/3*r;

        r = -1/h*(U_half(N+1,j)-U_half(N,j)+V_half(N,j+1)-V_half(N,j)) - Z(N,j);
        delta = r*h/3;
        V_half(N,j) = V_half(N,j) - delta;
        V_half(N,j+1) = V_half(N,j+1) + delta;
        U_half(N,j) = U_half(N,j) - delta;
        P(N,j) = P(N,j) + r;
        P(N,j+1) = P(N,j+1) - 1/3*r;
        P(N-1,j) = P(N-1,j) - 1/3*r;
        P(N,j-1) = P(N,j-1) - 1/3*r;
    end

    %四个角
    
    r = -1/h*(U_half(2,1)-U_half(1,1)+V_half(1,2)-V_half(1,1))-Z(1,1);
    delta = r*h/2;
    
    V_half(1,2) = V_half(1,2) + delta;
    U_half(2,1) = U_half(2,1) + delta;
    P(1,1) = P(1,1) + r;
    P(1,2) = P(1,2) - 1/2*r;
    P(2,1) = P(2,1) - 1/2*r;

    
    
    r = -1/h*(U_half(N+1,1)-U_half(N,1)+V_half(N,1+1)-V_half(N,1)) - Z(N,1);
    delta = r*h/2;
    
    V_half(N,1+1) = V_half(N,1+1) + delta;
    U_half(N,1) = U_half(N,1) - delta;
    P(N,1) = P(N,1) + r;
    P(N-1,1) = P(N-1,1) - 1/2*r;
    P(N,1+1) = P(N,1+1) - 1/2*r;
    
    
    r = -1/h*(U_half(2,N)-U_half(1,N)+V_half(1,N+1)-V_half(1,N)) - Z(1,N);
    delta = r*h/2;
    
    V_half(1,N) = V_half(1,N) - delta;
    U_half(2,N) = U_half(2,N) + delta;
    P(1,N) = P(1,N) + r;
    P(2,N) = P(2,N) - 1/2*r;
    P(1,N-1) = P(1,N-1) - 1/2*r;
    

    r = -1/h*(U_half(N+1,N)-U_half(N,N)+V_half(N,N+1)-V_half(N,N)) - Z(N,N);
    delta = r*h/2;
    
    V_half(N,N) = V_half(N,N) - delta;
    U_half(N,N) = U_half(N,N) - delta;
    P(N,N) = P(N,N) + r;
    P(N-1,N) = P(N-1,N) - 1/2*r;
    P(N,N-1) = P(N,N-1) - 1/2*r;


    



    
end




