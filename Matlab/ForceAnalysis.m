%% 受力分析，得运动微分方程
tic
%% 通用方程
%R---轨道半径
%r---不倒翁底部半径
%a0 a1 a2---摆角的位移，角速度，角加速度
%e---质心与圆心的距离
syms R r a0 a1 a2 e; 
thita = R*a0/(R-r); %底部半圆中心轴转过角度
phi = r*a0/(R-r);
aot = r*a2;  %圆心切向加速度
aon = r^2/(R-r)*(a1^2); %圆心法向加速度
aoct = e*a2; %质心相对圆心切向加速度
aocn= e*a1^2; %质心相对圆心法向加速度
acx = aot*cos(thita)-aon*sin(thita)-aoct; %质心x方向加速度
acy = aot*sin(thita)+aon*cos(thita)+aocn; %质心y方向加速度

%% 无电磁无激励
%Jc---不倒翁在质点C处惯性矩
%Ft---不倒翁与轨道接触点切向力
%Fn---不倒翁与轨道接触点法向力
%m---不倒翁质量
%g---重力加速度
syms Jc Ft Fn m g;
f11 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx == 0; %x方向力平衡
f21 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0; %y方向力平衡
eqns = [f11,f21]; 
vars = [Fn Ft];
[sol_Fn,sol_Ft] = solve(eqns,vars); %求解Fn，Ft
f31 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita); %力矩平衡


%% 无电磁向上激励
%Fe---激励力
syms ae;
f12 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx == 0;
f22 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0;
eqns = [f12,f22];
vars = [Fn Ft];
[sol_Fn,sol_Ft] = solve(eqns,vars);
f32 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)-Fe*((R-r)*sin(phi)-e*sin(a0));
f3new2 = collect(f32,Fe)-f31 %合并同类项
    
%% 有摩擦无电磁无激励
%Fs---摩擦力
%fs---滚动摩擦系数
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

%% 有电磁力无激励力
%Fm---电磁力
%Lc---电磁力与质心的距离
syms Fm Lc;
% f14 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-m*acx-2*Fm == 0; %x方向力平衡
% f24 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)-m*acy == 0; %y方向力平衡
% eqns = [f14,f24]; 
% vars = [Fn Ft];
% [sol_Fn,sol_Ft] = solve(eqns,vars); %求解Fn，Ft
% f34 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)+2*Fm*Lc; %力矩平衡
% f3new4 = collect(f34,Fm)-f31


%% 有电磁向上激励
% f1 = -Ft*cos(thita)-Fn*sin(thita)+m*g*sin(a0)-Fm-Fe*sin(a0)-m*acx == 0;
% f2 = -Ft*sin(thita)+Fn*cos(thita)-m*g*cos(a0)+Fe*cos(a0)-m*acy == 0;
% eqns = [f1,f2];
% vars = [Fn Ft];
% [sol_Fn,sol_Ft] = solve(eqns,vars)
% f3 = Jc*a2-sol_Ft*(r-e*cos(thita))+sol_Fn*e*sin(thita)+2*Fm*Lc-Fe*(r-e)*sin(a0);
% f3new4 = collect(f3,Fm) %合并同类项

%% 计算程序运行时间
toc