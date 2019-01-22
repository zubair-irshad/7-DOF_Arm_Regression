phi_mat=load('phifile_comp');

all_torque=load('dataTorque_sec.txt');

all_torque=all_torque((50+1):end,:);

phi_mat=phi_mat((50*7)+1:end,:);

[m,n]=size(phi_mat);

[o,p]=size(all_torque);

T=[];

for i = 1:o
    ta=all_torque(i,:)';
    tau=[T;ta];
    T=tau;
end


beta_param=load('betaparameters_hardware.txt'); 


% std_param is the numerical value of parameters with '1' being M1 '2'
% being 'MX1' and so on...

std_param=linspace(1,n,n);

r=rank(phi_mat);

%rr=phi_mat'*phi_mat;
%e=eig(rr);
%min_e=min(e);
%ridge=(phi_mat'*phi_mat)+(1e-5*eye(70));

%param=pinv(ridge)*phi_mat'*tau;

param=pinv(phi_mat'*phi_mat)*phi_mat'*tau;

theta_opt=param;

phi_mattest=load('phifile_comp_test');

all_torquetest=load('dataTorque_sec_test.txt');

[j,k]=size(all_torquetest);

time=linspace(1,j,j)';

tau_act=phi_mattest*theta_opt;

T1=zeros(j,k);

for i = 1:j
    ta1=tau_act(7*(i-1)+1:(7*(i-1)+7))';
    T1(i,:)=ta1;
end

tau_reshape=T1;


torque1=all_torquetest(:,1);
torque2=all_torquetest(:,2);
torque3=all_torquetest(:,3);
torque4=all_torquetest(:,4);
torque5=all_torquetest(:,5);
torque6=all_torquetest(:,6);
torque7=all_torquetest(:,7);

torque1_a=tau_reshape(:,1);
torque2_a=tau_reshape(:,2);
torque3_a=tau_reshape(:,3);
torque4_a=tau_reshape(:,4);
torque5_a=tau_reshape(:,5);
torque6_a=tau_reshape(:,6);
torque7_a=tau_reshape(:,7);

%diff_t=tau -tau_act;

diffparam=theta_opt-beta_param;

figure

subplot(2,3,1)
plot(time,torque1,'b',time,torque1_a,'g')
title('Joint1')

subplot(2,3,2)
plot(time,torque2,'b',time,torque2_a,'g')
title('Joint2')

subplot(2,3,3)
plot(time,torque3,'b',time,torque3_a,'g')
title('Joint3')

subplot(2,3,4)
plot(time,torque4,'b',time,torque4_a,'g')
title('Joint4')

subplot(2,3,5)
plot(time,torque5,'b',time,torque5_a,'g')
title('Joint5')

subplot(2,3,6)
plot(time,torque6,'b',time,torque6_a,'g')
title('Joint6')


