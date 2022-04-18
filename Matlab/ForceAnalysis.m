%% �������������˶�΢�ַ���
tic
%% ͨ�÷���
%R---����뾶
%r---�����̵ײ��뾶
%a0 a1 a2---�ڽǵ�λ�ƣ����ٶȣ��Ǽ��ٶ�
%e---������Բ�ĵľ���
syms R r a0 a1 a2 e; 
thita = R*a0/(R-r); %�ײ���Բ������ת���Ƕ�
phi = r*a0/(R-r);
aot = r*a2;  %Բ��������ٶ�
aon = r^2/(R-r)*(a1^2); %Բ�ķ�����ٶ�
aoct = e*a2; %�������Բ��������ٶ�
aocn= e*a1^2; %�������Բ�ķ�����ٶ�
acx = aot*cos(thita)-aon*sin(thita)-aoct; %����x������ٶ�
acy = aot*sin(thita)+aon*cos(thita)+aocn; %����y������ٶ�

%% �޵���޼���
%Jc---���������ʵ�C�����Ծ�
%Ft---�����������Ӵ���������
%Fn---�����������Ӵ��㷨����
%m---����������
%g---�������ٶ�
syms Jc Ft Fn m g;
f11 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx == 0; %x������ƽ��
f21 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0; %y������ƽ��
eqns = [f11,f21]; 
vars = [Fn Ft];
[sol_Fn,sol_Ft] = solve(eqns,vars); %���Fn��Ft
f31 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita); %����ƽ��


%% �޵�����ϼ���
%Fe---������
syms ae;
f12 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx == 0;
f22 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0;
eqns = [f12,f22];
vars = [Fn Ft];
[sol_Fn,sol_Ft] = solve(eqns,vars);
f32 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)-Fe*((R-r)*sin(phi)-e*sin(a0));
f3new2 = collect(f32,Fe)-f31 %�ϲ�ͬ����
    
%% ��Ħ���޵���޼���
%Fs---Ħ����
%fs---����Ħ��ϵ��
syms Fs fs;
% Fs = fs* Fn;
% f13 = -(Ft+Fs)*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx == 0;
% f23 = -(Ft+Fs)*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0;
% eqns = [f13,f23];
% vars = [Fn Ft];
% [sol_Fn,sol_Ft] = solve(eqns,vars);
% sol_Fs = fs*sol_Fn;
% f33 = Jc*a2-(sol_Ft+sol_Fs)*(r-e*cos(thita))+sol_Fn*e*sin(thita);
% f3new3 = collect(f33,fs)

%% �е�����޼�����
%Fm---�����
%Lc---����������ĵľ���
syms Fm Lc;
% f14 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx-2*Fm == 0; %x������ƽ��
% f24 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0; %y������ƽ��
% eqns = [f14,f24]; 
% vars = [Fn Ft];
% [sol_Fn,sol_Ft] = solve(eqns,vars); %���Fn��Ft
% f34 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)+2*Fm*Lc; %����ƽ��
% f3new4 = collect(f34,Fm)-f31


%% �е�����ϼ���
% f1 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-Fm-Fe*sin(a0)-m*acx == 0;
% f2 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)+Fe*cos(a0)-m*acy == 0;
% eqns = [f1,f2];
% vars = [Fn Ft];
% [sol_Fn,sol_Ft] = solve(eqns,vars)
% f3 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)+2*Fm*Lc-Fe*(r-e)*sin(a0);
% f3new4 = collect(f3,Fm) %�ϲ�ͬ����

%% �����������ʱ��
toc