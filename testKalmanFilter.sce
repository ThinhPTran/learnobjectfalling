// H = eye(2,2): This means that the measure readings are the same as the physical parameters we are considering.
// R =  eye(2,2)*MNV; this is the uncertainty of measurement. Constant in a short time.
// Q = eye(2,2)*PNV; this is the uncertainty of process. Constant in a short time.
// P = eye(2,2)*MNV: initial uncertainty of position and velocity is the uncertainty of measurement.
// Kalman filer example

clear all
clc

// Functions
function s = kalmanf(s)
    s.x = s.A*s.x + s.B*s.u;
    s.P = s.A *s.P * s.A' + s.Q;
    // Compute Kalman gain factor:
    K = s.P * s.H' * inv(s.H * s.P * s.H' + s.R);
    // Correction based on observation
    s.x = s.x + K*(s.z - s.H*s.x);
    s.P = s.P - K*s.H*s.P;
endfunction


// Main Program
// Constants
x0 = 0.0;
v0 = 0.0;
dt = 0.01;  
k = 10.0;
g = 10;
m = 1.0;
N = 200;

// Measurement noise variance
MNstd = 0.3;
MNV = MNstd*MNstd;

// Process noise variance
PNstd = 0.1;
PNV = PNstd*PNstd;

s = zeros(N,1)

t = [0:dt:(N-1)*dt]

// Dynamics of the system
s.A = [[1    dt];
           [0    (1.0-(k/m)*dt)]];
       
// Use control to include gravity (External controls)
s.B = eye(2,2); 
// Control matrix
s.u = [0    -g*m*dt]'; //Gravitational acceleration

// Initial state:
s.x = [x0 v0]';
s.P = eye(2,2)*MNV;
       
// Process noise covariance matrix
s.Q = eye(2,2)*PNV;

// Define measurement function to return the state
s.H = eye(2,2);

// Define a measurement error
s.R = eye(2,2)*MNV;

// Let's keep track of the noise by keeping detP
s.detP = det(s.P);  
s.z = zeros(2,1);

// Simulate falling in air, and watch the  filter track it
tru = zeros(N,2)
theo = zeros(N,2)
tru(1,:) = [x0 v0];
theo(1,:) = [x0 v0];
detP(1,:) = s.detP;

for iter = 2:N
    theo(iter,:) = (s(iter-1).A*theo(iter-1,:)'+s(iter-1).B*s(iter-1).u)';
    tru(iter,:) = (s(iter-1).A*tru(iter-1,:)'+s(iter-1).B*s(iter-1).u+PNstd*rand(2,1,'normal'))';
    s(iter-1).z = s(iter-1).H * tru(iter,:)' + MNstd*rand(2,1,'normal'); // create a meas.
    s(iter) = kalmanf(s(iter-1));  // perform a Kalman filter iteration
    detP(iter) = s(iter).detP; // keep track of "net" uncertainty
end

// Collect z
zvec = zeros(N,1);
xvec = zeros(N,1);
for iter = 1:N
    zvec(iter,1) = s(iter).z(1,1);
    xvec(iter,1) = s(iter).x(1,1); 
end

figure(1)
plot(t,tru(:,1),t,xvec,t,zvec,'*',t,theo(:,1))
legend('true','Kalman filtered','real','theory')



