clear all
clc


// Functions



// Main Program
// Constants
x0 = 0.0; 
v0 = 0.0;
dt = 0.01;
k = 10.0; 
g = 10; 
m = 1.0; 
N = 200; 

t = [0:dt:(N-1)*dt]
x = zeros(2,N)

x(1,1) = x0; 
x(2,1) = v0; 

u=[0 1]';

for iter = 2:N
    A=[[1    dt];
       [0    (1.0-(k/m)*dt)]];
    B=[[1   0];
       [0   -g*dt]]; 
    x(:,iter) = A*x(:,iter-1) + B*u;
end

figure(1)
plot(t,x(1,:),'r',t,x(2,:),'b')
legend('X','V');



