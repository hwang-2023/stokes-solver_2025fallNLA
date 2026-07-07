function [F, G, U_true, V_true] = initialize(N)

%计算初始的F,G与真解U_true,V_true
%其中, 把边界条件放在右端


h = 1/N;

F = zeros(N+1,N);

for j = 2:N-1
    for i = 2:N
        x = (i-1)*h;
        y = (j-1/2)*h;
        F(i,j) = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y) + x^2;
    end
end

for i = 2:N
    x = (i-1)*h;
    y = (1/2)*h;
    F(i,1) = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y) + x^2 + 1/h*(2*pi*(cos(2*pi*x)-1));

    x = (i-1)*h;
    y = (N-1/2)*h;
    F(i,N) = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y) + x^2 + 1/h*(-2*pi*(cos(2*pi*x)-1));
end


G = zeros(N,N+1);

for i = 2:N-1
    for j = 2: N
        x = (i-1/2)*h;
        y = (j-1)*h;
        G(i,j) = 4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x);
    end
end
for j = 2:N
    x = 1/2*h;
    y = (j-1)*h;
    G(1,j) = 4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x) + 1/h*(-2*pi*(cos(2*pi*y)-1));

    x = (N-1/2)*h;
    y = (j-1)*h;
    G(N,j) = 4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x) + 1/h*(2*pi*(cos(2*pi*y)-1));
end


U_true = zeros(N+1,N);


for i = 1:N+1
    for j = 1:N
        x = (i-1)*h;
        y = (j-1/2)*h;
        U_true(i,j) = (1-cos(2*pi*x))*sin(2*pi*y);
    end
end


V_true = zeros(N,N+1); 

for j = 1:N+1
    for i = 1:N
        x = (i-1/2)*h;
        y = (j-1)*h;
        V_true(i,j) = -(1-cos(2*pi*y))*sin(2*pi*x);
    end
end


end