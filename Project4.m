clc,clear
k = input('Enter 1 for Part A, 2 for Part B and 3 for Part C \n');
if k == 1
a = input(' Enter 1 For Calculating number of tube passes (n) \n Enter 2 for plotting the graph of number of tube passes (n) vs outlet temp of hot fluid \n')
if a == 1
T1 = 55; T2 = 15;
T3 = 2;
U = 700; n = zeros();
D = 0.015; L = 3;
Cp_h = 2630; Cp_c = 4186;
mh = 0.28; mc = 1.3;
A = n*pi*D*L;
Cmin = Cp_h*mh;
Cmax = Cp_c*mc;
c = Cmin/Cmax;
Qmax = Cmin*(T1-T3);
 Qact = Cp_h*mh*(T1-T2);
 eps = Qact/Qmax;
 P = (1+c^2)^0.5;
 R = (2/eps)-1-c;
 NTU = -1/P*log((R-P)/(R+P));
 n = NTU*Cmin/(U*pi*D*L);
 disp('Number of tube passes (n) comes out be');
 round(n) 
 Qact
end

if a == 2
T1 = 55; T2 = zeros(14,1);
T3 = 2; 
% deltaT1 = T1-T2;
R = zeros(14,1);
eps = zeros(14,1);
X = zeros(14,1);
Qact = zeros(14,1);
n = zeros(14,1);
NTU = zeros(14,1);
U = 700; 
D = 0.015; L = 3;
Cp_h = 2630; Cp_c = 4186;
mh = 0.28; mc = 1.3;
A = pi*D*L;
Cmin = Cp_h*mh;
Cmax = Cp_c*mc;
c = Cmin/Cmax;
 Qmax = Cmin*(T1-T3);
  P = (1+c^2)^0.5;
  n = 1;
 for i = 1:14
  NTU(i) = n*A*U/Cmin;
 X(i) = exp(-NTU(i).*P);
 R(i) = (1+X(i))./(1-X(i));
 eps(i) = 2/(1+c+P*R(i));
 Qact(i) = eps(i)*Qmax;
 T2(i) = T1 - Qact(i)./Cp_h;
 n = n+1;
 end
    n=1:14;
    figure
  plot(T2,n,'ob','LineWidth',3);
 hold on
 title('Graph of number of tube passes (n) vs outlet temp of hot fluid (Th,out)','FontSize',14);
 xlabel('Outlet temp of hot fluid (Th,out)','FontSize',13);
 ylabel(' Number of tube passes (n)','FontSize',13);
end

end

if k == 2 

T1 = 55; T2 = zeros(30,1);
T3 = 2; T4 = zeros(30,1);
R = zeros(30,1);
eps = zeros(30,1);
Qact = zeros(30,1);
n = zeros(30,1);
X = zeros(30,1);
NTU = zeros(30,1);
U = 700; 
D = 0.015; L = 3;
Cp_h = 2630; Cp_c = 4186;
mh = 0.28; mc = 1.3;
A = pi*D*L;
Cmin = Cp_h*mh;
Cmax = Cp_c*mc;
c = Cmin/Cmax;
 Qmax = Cmin*(T1-T3);
P = (1+c^2)^0.5;
n = 1;
for i = 1:30
 NTU(i) = n*A*U/Cmin;
 X(i) = exp(-NTU(i).*P);
 R(i) = (1+X(i))./(1-X(i));
 eps(i) = 2/(1+c+P*R(i));
 Qact(i) = eps(i)*Qmax;
 T2(i) = T1 - Qact(i)./Cp_h;
%  fprintf('Output temperature(Th,out) of hot  fluid(Ethanol) for the n = %i is = %f \n',n,T2(i));
 T4(i) = T3 + Qact(i)./Cp_c;
%  fprintf('Output temperature(Tc,out) of cold fluid(Water)   for the n = %i is = %f\n',n,T4(i));
%  fprintf('Effectiveness for the n = %i is = %f\n',n,eps(i));
n = n+1;
end

n = 1;
for i = 1:14
 fprintf('Output temperature(Th,out) of hot  fluid(Ethanol) for the n = %i is = %f \n',n,T2(i));
 fprintf('Output temperature(Tc,out) of cold fluid(Water)   for the n = %i is = %f\n',n,T4(i));
 fprintf('Effectiveness for the n = %i is = %f\n',n,eps(i));
 n = n+1;
end

figure
 n=1:30;
 plot(T4,n,'or','LineWidth',3);
hold on
title('Graph of number of tube passes (n) vs outlet temp of cold fluid (Tc,out)','FontSize',13);
xlabel('Outlet temp of cold fluid (Tc,out)','FontSize',14);
ylabel(' Number of tube passes (n)','FontSize',14);
end

if k == 3
    Cp_h = 2630; Cp_c = 4186;
mh = 0.28; mc = 1.3;
Cmin = Cp_h*mh;
Cmax = Cp_c*mc;
c = Cmin/Cmax;
U = 700;
D = 0.015; L = 3;
A = pi*D*L;

e = zeros(1,1);
NTU = U*A/Cmin;
a = input('Enter 1 for finding effectiveness of Parallel flow , Enter 2 for finding effectiveness of Counter flow ,Enter 3 for finding effectiveness of Cross flow(Cmax mixed) ,Enter 4 for finding effectiveness of Cross flow(Cmin mixed) \n');
if a == 1
e = (1-exp(-(NTU)*(1+c)))/(1+c);
e
end
if a == 2
e = (1-exp(-NTU*(1-c)))/(1-c*exp(-NTU*(1-c)));
e
end
if a == 3
  e = 1/c(1-exp(-c*(1-exp(-NTU))));
  e
end
if a == 4
   e = 1-exp(-(1/c)*(1-exp(-c*NTU)));
   e
end
end