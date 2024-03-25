clc
clear
close all
tic
T =100; N = 10^5; dt = T/N;
s = 0.5; k=120;
% 
% for k =[50,100,150,200]


%dB_k=sqrt(dt)*randn(N,k);
dB_k=sqrt(dt)*randn(N,k);
x0=0.69*ones(1,k); y0=0.01*ones(1,k);z0=0.1*ones(1,k);w0=0.2*ones(1,k);
%S0+I0+Q0+V0

% S0 = rand(1,k); I0 = rand(1,k); Q0 = rand(1,k);  V0 = 1-S0-I0-Q0;
% S0=1-randn(k,1)/4; I0=randn(k,1).*(1-S0)/4; Q0=randn(k,1).*(1-S0)/4;  V0=1-S0-I0-Q0;
% [S0 I0 Q0 V0]
% pause
Theta = 0;

A = 0.05; beTa = 0.7 ; theta_t = 0.05 ; gamma = 0.2; mu = 0.2 ; %epsilon_k =rand(1,k)*0.1+0.5
epsilon_k = rand(1,k)*0+0.2 ; q = 0.5 ; lambda = 0.4 ; h = 0.5 ;r=2.8;
phi = 0.001 ; alpha = 0.2
"'''''''p099"; D_k = rand(1,k)*0.5;   kt=[1:k];
%R0 = k_1n(2,n)*beTa/k_1n(1,n)*(gamma+mu+alpha)<1
%S0 = [0.8 ; 0.7] ; I0 = [0.1 ; 0.2]; Q0 = [0.05 ; 0.04] ; V0 = [0.05 ; 0.06] ;


x_em = zeros(N,k);
y_em = zeros(N,k);
z_em = zeros(N,k);
w_em = zeros(N,k);
SUM_em = zeros(N,k);

x_temp = x0; y_temp = y0; z_temp = z0; w_temp = w0;
for j=1:N
    x_kt = x_temp;
    y_kt = y_temp;
    z_kt = z_temp;
    w_kt = w_temp;
    Theta = sum(kt.*probabilityC(k,r).*(y_kt),'all')./k_1n(1,k,r); %theta_k       
    %Theta = Theta_i + Theta;
    x_temp = x_kt + (A-beTa*kt.*x_kt.*Theta-(A+q)*x_kt-epsilon_k.*x_kt+gamma*y_kt+alpha*x_kt.*y_kt+lambda*z_kt+h*w_kt)*dt - D_k.*x_kt.*y_kt.*dB_k(j,:);
    y_temp = y_kt + (beTa*kt.*x_kt.*Theta-(A+gamma+alpha)*y_kt+alpha*y_kt.*y_kt)*dt + D_k.*x_kt.*y_kt.*dB_k(j,:);
    z_temp = z_kt + (epsilon_k.*x_kt-(A+lambda)*z_kt+alpha*y_kt.*z_kt)*dt;%OK
    w_temp = w_kt + (q*x_kt-(A+h)*w_kt+alpha*y_kt.*w_kt)*dt;%OK

    digit_num = 15;

    x_temp = round(x_temp,digit_num);
    y_temp = round(y_temp,digit_num);
    z_temp = round(z_temp,digit_num);
    w_temp = round(w_temp,digit_num);

    temp = x_temp + y_temp + z_temp;
    w_temp = ones(1,k) - temp;
    
    x_em(j,:) = x_temp;
    y_em(j,:) = y_temp;
    z_em(j,:) = z_temp;
    w_em(j,:) = w_temp;

    SUM_em(j,:) = x_em(j,:) + y_em(j,:) + z_em(j,:) + w_em(j,:);
% 
%     if isempty(find(SUM_em(j,:)~=1)) == 0 % check SUM_em = 1
%         fprintf('Numerical error. The sum is not equal to 1 at T=%s.\n',num2str(j*dt));
%         t2 = fprintf('Please try to change the digit_num.\n');
%         break
%     end

%     AAA = alpha*y_kt.*(x_kt+y_kt+z_kt+w_kt);
%     BBB = alpha*y_kt;
%     if isequal(AAA,BBB)~=1
%         AAA-BBB
%         j
%         break
%     end
end

% p_y10 = plot( 0:dt:(T-100) , [y0(1,10);y_em(1:(length(0:dt:(T-100))-1),10)] ,'r'); hold on
% p_y50 = plot( 0:dt:(T-100) , [y0(1,50);y_em(1:(length(0:dt:(T-100))-1),50)] ,'b'); hold on
% p_y90 = plot( 0:dt:(T-100) , [y0(1,120);y_em(1:(length(0:dt:(T-100))-1),120)] ,'g'); hold on

p_y10 = plot( 0:dt:T , [y0(1,10);y_em(:,10)] ,'r'); hold on
p_y50 = plot( 0:dt:T , [y0(1,50);y_em(:,50)] ,'b'); hold on
p_y90 = plot( 0:dt:T , [y0(1,120);y_em(:,120)] ,'g'); hold on
%end
t = 0:dt:T;
R0 = (k_1n(2,k,r)*beTa)/(k_1n(1,k,r)*(A+gamma));
%RE=beTa/(gamma+alpha);

%fprintf('RE=%f\n',RE);
fprintf('R0=%f\n',R0);
It=10^(gamma+mu+alpha)*(R0-1);
log_I=(gamma+mu+alpha)*(R0-1);

xlabel('t');
ylabel('I_k');
%zlabel('w_k')
%axis([0,100,0,0.5]);
set(gca,'tickdir','out');
legend('k=10','k=50','k=120','location','south')
toc