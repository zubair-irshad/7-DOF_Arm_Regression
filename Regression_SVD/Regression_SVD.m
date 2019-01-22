
% % Load data (q,qdot,qdotdot,torque) available online and save it
% A=load('sarcos_inv.mat');
% 
% B=A.sarcos_inv;
% 
% all_q=B(1:1000,1:7);
% 
% all_dq=B(1:1000,8:14);
% 
% all_dqq=B(1:1000,15:21);
% 
% all_torque=B(1:1000,22:28);

phi_mat=load('phifile_comp');

all_torque=load('dataTorque_sec.txt');


all_torque=all_torque((50+1):end,:);

phi_mat=phi_mat((50*7)+1:end,:);

[m,n]=size(phi_mat);

[o,p]=size(all_torque);


beta_param=load('betaparameters_hardware.txt'); 

T=[];

for i = 1:o
    ta=all_torque(i,:)';
    tau=[T;ta];
    T=tau;
end

% 
% W=[];
% 
% T=[];
% 
% %Load data from simulation
% 
% all_q=load('dataQ_train.txt');
% all_dq=load('dataQdot_train.txt');
% all_dqq=load('dataQdotdot_train.txt');
% all_torque=load('dataTorque_train.txt');
% 
% %Unwrap q,dq,dqq,all_torque into understandable parameters
% 
% for i = 1:2000
%    q=all_q(i,:);
%    dq=all_dq(i,:);
%    dqq=all_dqq(i,:);
%    
%    ph=PHI(q,dq,dqq);
%    phi_mat=[W;ph];
%    
%    W=phi_mat;
%    
%    
%    ta=all_torque(i,:)';
%    
%    tau=[T;ta];
%    
%    T=tau;
%    
% end
% 
% [m,n]=size(phi_mat);
% 
% filename = 'param_urdf.xlsx';
% Z = xlsread(filename);
% Z1=reshape(Z,[11,7]);
% 
% 
% 
% 
% for i = 1:7
%     s=[Z1(2,i), Z1(3,i), Z1(4,i)];
%     iMat=[Z1(5,i), Z1(6,i), Z1(7,i),Z1(8,i), Z1(9,i), Z1(10,i)];
%     tMat=[(s(2)*s(2)+s(3)*s(3)), (-s(1)*s(2)),(-s(1)*s(3)),(s(1)*s(1)+s(3)*s(3)),(-s(2)*s(3)),(s(1)*s(1)+s(2)*s(2))];
%     pm1=Z1(1,i);
%     
%     iMat1= iMat +pm1*tMat;
%     
%     
%     Z1(5:10,i)=iMat1;
%     
% end
% 
% 
% 
% Z2=reshape(Z1,[],1);  
% 
% tau_ref=phi_mat*Z2;


% std_param is the numerical value of parameters with '1' being M1 '2'
% being 'MX1' and so on...



std_param=linspace(1,n,n);

r=rank(phi_mat);



[U,S,V]=svd(phi_mat);



for i =1:n
    tao=1000*eps*S(1,1);
    if S(i,i)<=tao
        S(i,i)=0;
    end    
end

rank_S=rank(S);

tau_ref=phi_mat*beta_param;


V1= V(:,1:rank_S);
S1=S(1:rank_S,1:rank_S);
U1= U(:,1:rank_S);


theta_opt=V1*inv(S1)*U1'*tau;

theta_opt_fin = beta_param + (V1*inv(S1)*U1'*(tau-tau_ref));

diff_beta=theta_opt - beta_param;





% Testing set
%Uncoment for testing set%
% 
% 
% all_q_test=load('dataQ_test.txt');
% all_dq_test=load('dataQdot_test.txt');
% all_dqq_test=load('dataQdotdot_test.txt');
% all_torque_test=load('dataTorque_test.txt');
% 
% 
% time=load('dataTime_test.txt');

%Plot simulation torques



% W_test=[];
% 
% T_test=[];
% 
% [s1,s2]=size(all_q_test);
% 
% for i = 1:s1
%    q_test=all_q_test(i,:);
%    dq_test=all_dq_test(i,:);
%    dqq_test=all_dqq_test(i,:);
%    
%    ph_test=PHI(q_test,dq_test,dqq_test);
%    phi_mat_test=[W_test;ph_test];
%    
%    W_test=phi_mat_test;
%    
%    
%    ta_test=all_torque_test(i,:)';
%    
%    tau_test=[T;ta_test];
%    
%    T_test=tau_test;
%    
% end





%tau_act=phi_mat_test*theta_opt;

%tau_reshape=reshape(tau_act,[s1,s2]);


% check it on testing set


phi_mattest=load('phitest.txt');

all_torquetest=load('dataTorque_RHStest.txt');

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


%%%

% [Q,R]=qr(phi_mat);
% 
% k=diag(R);
% 
% indices=zeros(n,1);
% 
% % count number of zeros in indices (independent parameters)
% 
% nindependant=sum(indices==0);
% 
% ndependant=size(k)-size(indices);
% 
% for i =1:size(k)
%     tao=1000*eps*abs(max(k));
%     if abs(k(i))<=tao
%         indices(i)=1;
%     end    
% end

%%%

% ic=zeros(n,1);
% 
% for i =1:n
%     if norm(phi_mat(:,i)) == 0;
%         ic(i)=1;
%     end
% end
% 
% un_ident=find(ic==1);
% 
% phi_mat1=phi_mat;
% 
% phi_mat(:,un_ident)=[];
% 
% std_param(:,un_ident)=[];
% 
% 
% 
% ident=[];
% ind_dep=[];
% 
% [o,p]=size(phi_mat);
% 
% L=zeros(o,p);
% 
% phi_matprior=phi_mat;
% 
% for i = 1:p
%     
%     fci=phi_mat(:,i);
%     
%     phi_mat(:,i)=[];
%     
%     rwi=rank(phi_mat);
%     
%     if rwi <r
%         ident(end+1)=std_param(i);
%         L(:,i)=fci;
%         
%     else
%         ind_dep(end+1)=i;
%         
%     end
%         
%     if i==1
%         
%         B=[fci phi_mat];
%         
%         phi_mat=B;
%         
%     end
%     
%     if i>1 && i<n
%         B=[phi_mat(:,1:i-1) fci phi_mat(:,i:p-1)];
%         
%         phi_mat=B;
%         
%     end
%     
%     if i==n
%         B=[phi_mat fci];
%         phi_mat=B;
%     end
%     
% end
% 
% %%%%
% 
% % ic2=zeros(p,1);
% % 
% % for i =1:p
% %     if norm(L(:,i)) == 0;
% %         ic2(i)=1;
% %     end
% % end        
% % 
% % ind=find(ic2==1);
% 
% 
% %%%%
% 
% L(:,ind_dep)=[];
% 
% DBS=ident;
% 
% rank_L=rank(L);
% 
% [~,q]=size(ind_dep);
% 
% ident_LC=[];
% 
% 
% for i =1:q
%     
%     rank_L=rank(L);
%     
%     fc=phi_mat(:,ind_dep(i));
%     
%     L1=[L fc];
%     
%     rank_L1=rank(L1);
%     
%     if rank_L1>rank_L
%         
%         DBS(end+1)= std_param(ind_dep(i));
%    
%         L(:,end+1)=fc;
%     
%     else
%         
%         ident_LC(end+1)=std_param(ind_dep(i));
%         
%     end
%     
%     
% end
% 
% s_LC=size(ident_LC);
% LC_theta=zeros(s_LC);
% 
% ident_theta=pinv(L'*L)*L'*tau;
% 
% for i = 1:s_LC
%     
%     fc_LC=phi_mat1(:,ident_LC(i));
%     
%     alpha=pinv(L'*L)*L'*fc_LC;
%     
%     LC_theta(i)=ident_theta'*alpha; 
% end
% 
% i_theta=[ident_theta DBS'];
% li_theta=[LC_theta' ident_LC']; 
% 
% 
% 
% % alpha = 0.01;
% % num_iters = 400;
% % 
% % theta = zeros(77, 1);
% % [theta, J_history] = gradientDescentMulti(phi_base, tau, theta, alpha, num_iters);
% % 
