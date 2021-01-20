clc; clear; warning('off');
f=1046.502; %cutoff frequency
w=2*pi*f;  %cutoff freq in radians
dt= 29.89e-6; %sampling time
s =sym(1);
omega = linspace(0,2, 1001); %defining omega
z = exp(1i*omega); %defining Z

%%
%Setting up our poles
s1=w*exp(1i*5*pi/8);
s4=w*exp(1i*11*pi/8);
s2=w*exp(1i*7*pi/8);
s3=w*exp(1i*9*pi/8);

% Part1.A 
Z = [];
P = [-2.5163e+03+6.0748e+03i -2.5163e+03-6.0748e+03i -6.0748e+03+2.5163e+03i -6.0748e+03-2.5163e+03i];
K = abs(5*s1*s4*s3*s2);
Hcs = zpk(Z,P,K);
tf(Hcs);
Hc_s = tf(Hcs); %obtaining our correct Hc_s
DCGAIN = dcgain(Hc_s);%verifying DC gain


%%
% Part1.B & Part1.C
% discretization of the system
Zexact=log(z)/dt; %our Z-domain exact equation
Ztustin=(2/dt)*((z-1)./(z+1));% our Z-domain Tustin Approximation equation

HdExactE=(9.345e15)./((Zexact+(-2.5163e+03+6.0748e+03i)).*(Zexact+(-2.5163e+03-6.0748e+03i)).*(Zexact+(-6.0748e+03+2.5163e+03i)).*(Zexact+(-6.0748e+03-2.5163e+03i))); %Hd_exact(z)
HdTustin=(9.345e15)./((Ztustin+(-2.5163e+03+6.0748e+03i)).*(Ztustin+(-2.5163e+03-6.0748e+03i)).*(Ztustin+(-6.0748e+03+2.5163e+03i)).*(Ztustin+(-6.0748e+03-2.5163e+03i)));%Hd_tustin(z)
% semilogx(omega/(2*pi*dt),(abs(HdTustin-HdExactE))) %determined our absolute error
%%
% Part1.E

% semilogx(omega/dt, 20*log10(abs(HdExactE))) %plot Hd_exact(z)
% semilogx(omega/dt, 20*log10(abs(HdTustin))) %plot Hd_tustin(z)

%figure(1)
%semilogx(omega/dt,20*log10(abs(HdTustin)))
%title('H_d(Z)')
%xlim([0 8000])
%ylabel('dB')
%xlabel('Frequency rad/s')
%%
% Part2.B
% Determine a Parallel-Form Implementation
Hd_tustin = c2d(Hc_s, dt, 'tustin');
[num_HdT, den_HdT] = tfdata(Hd_tustin);
num_HdT= num_HdT{1};
den_HdT= den_HdT{1};

[A_mat,B_mat,C_mat,D_mat] = tf2ss(num_HdT,den_HdT); %matlab does a weird thing by being upside down or something, so I hardcoded the values to the correct form

%A= [0 1 0 0; 0 0 1 0; 0 0 0 1; 3.4880 -4.5920 2.7010 -.5988];
%B= [0; 0; 0; 1];
%C= [0.0001 0.0024 0.0005 0.0027];
%D= 3.6080e-04;

[x1,y1]= ss2tf(A_mat, B_mat, C_mat, D_mat);
system3= tf(x1,y1,dt); % Parallel Form SS Rep -> Tranfer Function
%figure(1)
%pzmap(system3)

[v, d]= eig(A_mat);%Getting the Eigen Vectors of the A matrix..

P = [real(v(:,1)), imag(v(:,1)), real(v(:,3)) imag(v(:,3))]; % Creating our P matrix

%Multiplying our ABC by the P matrix
A_parallel =  inv(P)*A_mat * P;
B_parallel =  inv(P)*B_mat;
C_parallel =  C_mat * P;

A1 = A_parallel(1:2,1:2);
B1 = B_parallel(1:2);
C1 = C_parallel(1:2);
[n1, d1]= ss2tf(A1, B1, C1, 0);

A2 = A_parallel(3:4,3:4);
B2 = B_parallel(3:4);
C2 = C_parallel(3:4);
[n2, d2]= ss2tf(A2, B2, C2, 0);
%%
%Part2.C
%Approx the multiplierse of your filter with your own fix point rep.

%PZMap of our original Hd_tustin
%figure(2)
[g1,h1]= ss2tf(A_parallel, B_parallel, C_parallel, 0);
system1= tf(g1,h1,dt);
%pzmap(Hd_tustin)

%PZMap of our Parallel-Form Space State Representation
%figure(3)
%pzmap(system1)

%[n -f] = [32, -16]
A_FixPoint = fi(A_parallel, 1, 32, 16);
B_FixPoint = fi(B_parallel, 1, 32, 16);
C_FixPoint = fi(C_parallel, 1, 32, 16);

A_Fix = [ 0.9129 0.1674 0 0; -0.1674 0.9129 0 0; 0 0 0.8314 0.0631; 0 0 -0.0631 0.8314];
B_Fix = [ -280.4633; -178.6837; 246.9775; 735.9525];
C_Fix = [ 0.0026 0.0009 0.0025 0.0004];

%Poles and Zeros of the Fix Point Approximation
%figure(4)
[g2,h2] = ss2tf(A_Fix,B_Fix,C_Fix,0);
system2 = tf(g2,h2,dt);
%pzmap(system2);

%Finding the Distance Between the Poles
Parallel_Poles = pole(system1);

Approx_Poles = pole(system2);
Difference_Poles= Parallel_Poles - Approx_Poles;



%%
%Part2.D
% Frequency Response

%Sys2FreqResponse= (.0218*z.^3 - 0.05613 *z.^2 + 0.05827*z -0.01823)./(z.^4 -3.489 * z.^3 + 4.593 * z.^2 - 2.702 *z + 0.5989);
%figure(1)
%semilogx(omega/dt, 20*log10(abs(Sys2FreqResponse)))
%title('Frequency Response of Fix Point')
%ylabel('dB')
%xlabel('Frequency rad/s')

%% Part 3
% Part 3.a

% We did not obtain any feedback, therefore, we continued as if part 2 was
% fully correct.

%%
% Part 3.B

%A_parallel Matrix
A = [ 0.9129    0.1674         0         0;
    -0.1674     0.9129         0         0;
           0          0    0.8314    0.0631;
           0          0   -0.0631    0.8314];
%B_Parallel Matrix     
B = [ -280.4634;
     -178.6387;
      246.9775;
      735.9525];
%C_Parallel Matrix  
C = [0.0026  0.0009  0.0025  0.0004];
  
% Identity Matrix
B1 = [1; 0; 0; 0];
B2 = [0; 1; 0; 0];
B3 = [0; 0; 1; 0];
B4 = [0; 0; 0; 1];

[num,den]=ss2tf(A,B,C,0);
[num1,den1]=ss2tf(A,B1,C,0);
[num2,den2]=ss2tf(A,B2,C,0);
[num3,den3]=ss2tf(A,B3,C,0);
[num4,den4]=ss2tf(A,B4,C,0);

%Partial Fraction Expansion
[r1,p1,k1] = residue(num,den);
[r2,p2,k2] = residue(num1,den1);
[r3,p3,k3] = residue(num2,den2);
[r4,p4,k4] = residue(num3,den3);
[r5,p5,k5] = residue(num4,den4);

%Cartesian to Polar 
[rho_r1, theta_r1] = cart2pol(real(r1), imag(r1));
[rho_r2, theta_r2] = cart2pol(real(r2), imag(r2));
[rho_r3, theta_r3] = cart2pol(real(r3), imag(r3));
[rho_r4, theta_r4] = cart2pol(real(r4), imag(r4));
[rho_r5, theta_r5] = cart2pol(real(r5), imag(r5));
[rho_p1, theta_p1] = cart2pol(real(p1), imag(p1));
[rho_p2, theta_p2] = cart2pol(real(p2), imag(p2));
[rho_p3, theta_p3] = cart2pol(real(p3), imag(p3));
[rho_p4, theta_p4] = cart2pol(real(p4), imag(p4));
[rho_p5, theta_p5] = cart2pol(real(p5), imag(p5));

%Lnorms
L1 = abs(2.*rho_r1(1)./rho_p1(1) + 2.*rho_r1(3)./(1 - rho_p1(3)));
L2 = abs(2.*rho_r2(1)./rho_p2(1) + 2.*rho_r2(3)./(1 - rho_p2(3)));
L3 = abs(2.*rho_r3(1)./rho_p3(1) + 2.*rho_r3(3)./(1 - rho_p3(3)));
L4 = abs(2.*rho_r4(1)./rho_p4(1) + 2.*rho_r4(3)./(1 - rho_p4(3)));
L5 = abs(2.*rho_r5(1)./rho_p5(1) + 2.*rho_r5(3)./(1 - rho_p5(3)));

%Bounds
M = (2^-7)/6;

%Fraction Bits
fu = -log2(M/L1);
fy = -log2(M/4);
f1 = -log2(M/(3*L2));
f2 = -log2(M/(3*L3));
f3 = -log2(M/(3*L4));
f4 = -log2(M/(3*L5));

%Fraction Bit Outputs
fprintf("Fraction bits for  input: %4.0f\n", -log2(M/L1));
fprintf("Fraction bits for  output: %4.0f\n", -log2(M/4));
fprintf("Fraction bits for  F1: %4.0f\n", -log2(M/(3*L2)));
fprintf("Fraction bits for  F2: %4.0f\n", -log2(M/(3*L3)));
fprintf("Fraction bits for  F3: %4.0f\n", -log2(M/(3*L4)));
fprintf("Fraction bits for  F4: %4.0f\n", -log2(M/(3*L5)));

%Part 3.C
    %See Report
    
%Part 3.D

%[18,-16] representation for Input
Linput = 2;

%Gu -> input, x1, x2, x3, x3 to 
guu = abs(2*(rho_r1(1))/(1 - rho_p1(1)) + 2*(rho_r1(3))/(rho_p1(3)));
gux1 = abs(2*(rho_r2(1))/(1 - rho_p2(1)) + 2*(rho_r1(3))/(rho_p1(3)));
gux2 = abs(2*(rho_r3(1))/(1 - rho_p3(1)) + 2*(rho_r1(3))/(rho_p1(3)));
gux3 = abs(2*(rho_r4(1))/(1 - rho_p4(1)) + 2*(rho_r1(3))/(rho_p1(3)));
gux4 = abs(2*(rho_r5(1))/(1 - rho_p5(1)) + 2*(rho_r1(3))/(rho_p1(3)));


ir1b = abs(B(1) * Linput);
ir2b = abs(B(2) * Linput);
ir3b = abs(B(3) * Linput);
ir4b = abs(B(4) * Linput);

igx1 = abs(gux1 * Linput);
ir11 = abs(A(1,1) * igx1);
ir21 = abs(A(2,1) * igx1);
ir1c = abs(C(1) * igx1);

igx2 = abs(gux2 * Linput);
ir12 = abs(A(1,2) * igx2);
ir22 = abs(A(2,2) * igx2);
ir2c = abs(C(2) * igx2);

igx3 = abs(gux3 * Linput);
ir43 = abs(A(4,3) * igx3);
ir33 = abs(A(3,3) * igx3);
ir3c = abs(C(3) * igx3);

igx4 = abs(gux4 * Linput);
ir34 = abs(A(3,4) * igx4);
ir44 = abs(A(4,4) * igx4);
ir4c = abs(C(4) * igx4);

[num1s,den1s]=ss2tf(A,B,[C(1) 0 0 0],B(1));
[num2s,den2s]=ss2tf(A,B,[0 C(2) 0 0],B(2));
[num3s,den3s]=ss2tf(A,B,[0 0 C(3) 0],B(3));
[num4s,den4s]=ss2tf(A,B,[0 0 0 C(4)],B(4));

[r1s, p1s, k1s] = residue(num1s,den1s);
[r2s, p2s, k2s] = residue(num1s,den1s);
[r3s, p3s, k3s] = residue(num1s,den1s);
[r4s, p4s, k4s] = residue(num1s,den1s);

[rho_r1s, theta_r1s] = cart2pol(real(r1s), imag(r1s));
[rho_r2s, theta_r2s] = cart2pol(real(r2s), imag(r2s));
[rho_r3s, theta_r3s] = cart2pol(real(r3s), imag(r3s));
[rho_r4s, theta_r4s] = cart2pol(real(r4s), imag(r4s));

[rho_p1s, theta_p1s] = cart2pol(real(p1s), imag(p1s));
[rho_p2s, theta_p2s] = cart2pol(real(p2s), imag(p2s));
[rho_p3s, theta_p3s] = cart2pol(real(p3s), imag(p3s));
[rho_p4s, theta_p4s] = cart2pol(real(p4s), imag(p4s));

ir1s = abs((B(1) + ((2*rho_r1s(1))/(1 - rho_p1s(1)))) * Linput);
ir2s = abs((B(2) + ((2*rho_r2s(2))/(1 - rho_p2s(2)))) * Linput);
ir3s = abs((B(3) + ((2*rho_r3s(3))/(1 - rho_p3s(3)))) * Linput);
ir4s = abs((B(4) + ((2*rho_r4s(4))/(1 - rho_p3s(3)))) * Linput);

irp1 = abs(ir1c + ir2c);
irp3 = abs(ir3c + ir4c);

[numr1,denr1]=ss2tf(A,B,[C(1) 0 0 0],B(1));
[numr2,denr2]=ss2tf(A,B,[0 C(2) 0 0],B(2));
[numr3,denr3]=ss2tf(A,B,[0 0 C(3) 0],B(3));
[numr4,denr4]=ss2tf(A,B,[0 0 0 C(4)],B(4));

%Integer Bit Outputs
fprintf("int bits for  Linput: %4.0f\n", log2(Linput) + 2)
fprintf("int bits for  r1b: %4.0f\n", log2(ir1b) + 2)
fprintf("int bits for  r2b: %4.0f\n", log2(ir2b) + 2)
fprintf("int bits for  r3b: %4.0f\n", log2(ir3b) + 2)
fprintf("int bits for  r4b: %4.0f\n", log2(ir4b) + 2)
fprintf("int bits for  gx1: %4.0f\n", log2(igx1) + 2)
fprintf("int bits for  r11: %4.0f\n", log2(ir11) + 2)
fprintf("int bits for  r21: %4.0f\n", log2(ir21) + 2)
fprintf("int bits for  r1c: %4.0f\n", log2(ceil(ir1c)) + 2)
fprintf("int bits for  gx2: %4.0f\n", log2(igx2) + 2)
fprintf("int bits for  r12: %4.0f\n", log2(ir12) + 2)
fprintf("int bits for  r22: %4.0f\n", log2(ir22) + 2)
fprintf("int bits for  r2c: %4.0f\n", log2(ceil(ir2c)) + 2)
fprintf("int bits for  gx3: %4.0f\n", log2(igx3) + 2)
fprintf("int bits for  r43: %4.0f\n", log2(ir43) + 2)
fprintf("int bits for  r33: %4.0f\n", log2(ir33) + 2)
fprintf("int bits for  r3c: %4.0f\n", log2(ceil(ir3c)) + 2)
fprintf("int bits for  gx4: %4.0f\n", log2(igx4) + 2)
fprintf("int bits for  r34: %4.0f\n", log2(ir34) + 2)
fprintf("int bits for  r44: %4.0f\n", log2(ir44) + 2)
fprintf("int bits for  r4c: %4.0f\n", log2(ceil(ir4c)) + 2)
fprintf("int bits for  r1s: %4.0f\n", log2(ir1s) + 2)
fprintf("int bits for  r2s: %4.0f\n", log2(ir2s) + 2)
fprintf("int bits for  r3s: %4.0f\n", log2(ir3s) + 2)
fprintf("int bits for  r4s: %4.0f\n", log2(ir4s) + 2)
fprintf("int bits for  rp1: %4.0f\n", log2(ceil(irp1)) + 2)
fprintf("int bits for  rp3: %4.0f\n", log2(ceil(irp3)) + 2)
fprintf("int bits for  r1: %4.0f\n", log2(igx1) + 2)
fprintf("int bits for  r2: %4.0f\n", log2(igx2) + 2)
fprintf("int bits for  r3: %4.0f\n", log2(igx3) + 2)
fprintf("int bits for  r4: %4.0f\n", log2(igx4) + 2)















