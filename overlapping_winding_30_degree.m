% clear
% clc
%% Institude: HuaZhong University of Science and Technology
% 机构： 华中科技大学电气与电子工程学院
%% Written by Xu Shouyu
% 作者： 徐首彧
% 指导老师： 叶才勇
%% 输入电机尺寸及参数

Rsb=0.16;                  % 定子槽底部半径
Rt=0.115;                  % 定子槽顶部半径/定子齿部半径
Rs=0.113;                  % 定子内半径
Rm=0.1025;                 % 磁极半径
Rr=0.0895;                 % 转子半径
boa=0.02655;               % boa是槽开口角
bsa=0.0607;                % bsa是槽宽角
ap=0.833;                  % 极弧系数
Rsm=(Rt+Rsb)/2;            % 槽内平均半径

p=2;                     % 极对数
Ns=48;                   % 齿槽个数
theta_s=2*pi/Ns;         % 定子齿槽极距，两个齿槽之间的间距，与选取半径无关，有槽数有关
                         % 第一个定子槽中心角的位置设为theta_s/2

a0=0;

u0=4*pi*1e-7;  % 真空磁导率
ur=1.05;       % PM的相对磁导率，解析模型的假设
Br=0.4;          % PM的剩磁/T，有限元中使用钕铁硼材料的剩磁密度

rpm=10000;     % rpm为电机转速
f=rpm*p/60;    % f为电枢电流频率
T=1/f;         % T为转子转一周的周期
wr=2*pi*f;     % wr为电角频率
we=wr/p;       % we为机械角频率
t=0;           % t为时刻

N=100;         % 气隙及磁极谐波次
M=100;         % 齿槽谐波次数取值数取值。在这组参数下，最大谐波次数取值。减小Rr,Rm,Rs半径间的差别，影响不大
ii=100;
K=(2*ii-1)*p;

%% 选择充磁方式 1.径向充磁 2.平行充磁

Mrck=zeros(K,1);
Mrsk=zeros(K,1);
Mack=zeros(K,1);
Mask=zeros(K,1);
magnetization_method=1;

switch(magnetization_method)
    case 1
        % radial magnetization 径向磁化
        for i = 1:ii
            n=2*i-1;
            k=n*p;

            Mrk=(4*p*Br/(k*pi*u0))*sin(k*pi*ap/(2*p));
            Mak=0;

            Mrck(k,1)=Mrk*cos(k*we*t+k*a0);
            Mrsk(k,1)=Mrk*sin(k*we*t+k*a0);
            Mack(k,1)=-Mak*sin(k*we*t+k*a0);
            Mask(k,1)=Mak*cos(k*we*t+k*a0);
        end

    case 2
        % parallel magnetization 平行磁化 k=1时A2k的分母为0要注意 mnjui8
        for i = 1:ii
            n=2*i-1;
            k=n*p;

            A1k=sin((k+1)*ap*(pi/(2*p)))/((k+1)*ap*(pi/(2*p)));
            if k == 1
                A2k=1;
            else
                A2k=sin((k-1)*ap*(pi/(2*p)))/((k-1)*ap*(pi/(2*p)));
            end

            Mrk=(Br/u0)*ap*(A1k+A2k);
            Mak=(Br/u0)*ap*(A1k-A2k);

            Mrck(k,1)=Mrk*cos(k*we*t+k*a0);
            Mrsk(k,1)=Mrk*sin(k*we*t+k*a0);
            Mack(k,1)=-Mak*sin(k*we*t+k*a0);
            Mask(k,1)=Mak*cos(k*we*t+k*a0);
        end
end


%% 选择绕组在槽内的空间分布 
% 移相角为30 degree

I_phase_max=693;                                   % 相电流幅值
m_phase=12;                                        % 绕组相数，这里绕组相数为12
Lys=2;                                             % 绕组层数，这里绕组层数为2
Nw=(Ns*Lys)/(2*m_phase);                           % 每相绕组的总线圈数，这里一相绕组有4个线圈
theta0=(178/180)*pi;
% S_slot=4.17e-4/2;
S_slot=93.57e-6;
t_region=linspace(0,T,1000);

% 为十二相绕组的相电流赋值
% I_phase(1,1)，I_phase(5,1)，I_phase(9,1) 分别是A1,B1,C1相绕组
% I_phase(2,1)，I_phase(6,1)，I_phase(10,1)分别是A2,B2,C2相绕组
% I_phase(3,1)，I_phase(7,1)，I_phase(11,1)分别是A3,B3,C3相绕组
% I_phase(4,1)，I_phase(8,1)，I_phase(12,1)分别是A4,B4,C4相绕组
I_phase=zeros(1000,m_phase);                              % 存储十二相电机的每相电流
for num_phase= 1:12
    for ti = 1:1000
        t=t_region(1,ti);
        I_phase(ti,num_phase)=I_phase_max*sin(wr*t+theta0-(num_phase-1)*pi/6);
    end
end

% 绘制十二相电流波形
figure
for i = 1:m_phase
    if i >= 1 && i <= 4
        NAME=["A",i,"相"];
        plot(t_region,I_phase(:,i),'-r','DisplayName',char(NAME))
        hold on
    elseif i >= 5 && i <= 8
        NAME=["B",i-4,"相"];
        plot(t_region,I_phase(:,i),'--b','DisplayName',char(NAME))
        hold on
    elseif i >= 9 && i <= 12
        NAME=["C",i-8,"相"];
        plot(t_region,I_phase(:,i),'-.g','DisplayName',char(NAME))
        hold on
    end
end

title("十二相电流波形")
xlabel("时间 / s")
ylabel("电流 / A")
set(gca,'FontSize',14)
legend

% J_phase为十二相绕组在槽内的面电流密度，J_low为下层绕组的排列，J_up为上层绕组的排列
I_phase=I_phase(1,:)';
% J_up=zeros(1000,48);
% J_low=zeros(1000,48);
% 
% figure
% for m = 1:m_phase
%     for nw = 1:Nw
%         % 假如相数为[1,2,5,6,9,10]，也就是A1B1C1，A2B2C2相，
%         if ismember(m,[1,2,5,6,9,10])
%             up_num=1+2*(m-1)+12*(nw-1);
%             low_num=up_num+11;
% 
%             if up_num==48
%             elseif up_num>48
%                 up_num=up_num-48;
%             end
%             if low_num==48
%             elseif low_num>48
%                 low_num=low_num-48;
%             end
% 
%             ai_up=theta_s/2+(up_num-1)*pi/24;                            % 上层绕组的角度
%             ai_low=theta_s/2+(low_num-1)*pi/24;                          % 下层绕组的角度
% 
%             % 相A1B1C1,A2B2C2的符号是红色的
%             if (-1)^(nw+1) == 1
%                 up_sign='or';
%                 low_sign='xr';
%             else
%                 up_sign='xr';
%                 low_sign='or';
%             end
% 
%             plot((Rsm+Rsb)/2*cos(ai_up),(Rsm+Rsb)/2*sin(ai_up),up_sign);  % 画上层绕组,半径为(Rsm+Rsb)/2
%             hold on
%             plot((Rt+Rsm)/2*cos(ai_low),(Rt+Rsm)/2*sin(ai_low),low_sign); % 画下层绕组,半径为(Rt+Rsm)/2
%             hold on
% 
%         %假如相数为[3,4,7,8,11,12]
%         elseif ismember(m,[3,4,7,8,11,12])
%             low_num=1+2*(m-1)+12*(nw-1);
%             up_num=low_num+11;
% 
%             if up_num==48
%             elseif up_num>48
%                 up_num=up_num-48;
%             end
%             if low_num==48
%             elseif low_num>48
%                 low_num=low_num-48;
%             end
% 
%             ai_up=theta_s/2+(up_num-1)*pi/24;     % 上层绕组的角度
%             ai_low=theta_s/2+(low_num-1)*pi/24;   % 下层绕组的角度
% 
%             % 相A3B3C3,A4B4C4的符号是蓝色的
%             if (-1)^(nw+1) == 1
%                 up_sign='xb';
%                 low_sign='ob';
%             else
%                 up_sign='ob';
%                 low_sign='xb';
%             end
% 
%             plot((Rsm+Rsb)/2*cos(ai_up),(Rsm+Rsb)/2*sin(ai_up),up_sign);  % 画上层绕组,半径为(Rsm+Rsb)/2
%             hold on
%             plot((Rt+Rsm)/2*cos(ai_low),(Rt+Rsm)/2*sin(ai_low),low_sign); % 画下层绕组,半径为(Rt+Rsm)/2
%             hold on
% 
%         end
% 
%         J_up(:,up_num)=J_phase(:,m)*((-1)^(nw+1));
%         J_low(:,low_num)=J_phase(:,m)*((-1)^nw);
% 
%     end
% end
% xlim([-0.2 0.2]);
% ylim([-0.2 0.2]);
% axis square
% 
% Jt1=J_low(1,1:48)';
% Jt2=J_up(1,1:48)';
% 
% It1=Jt1*(S_slot/2);
% It2=Jt2*(S_slot/2);

S_slot_up=(bsa/2)*(Rsb^2-Rsm^2);     % 上层槽的面积
S_slot_down=(bsa/2)*(Rsm^2-Rt^2);    % 下层槽的面积

It1=[-I_phase(7,1) ; I_phase(2,1) ; -I_phase(8,1) ; -I_phase(9,1) ; I_phase(3,1) ; -I_phase(10,1) ; I_phase(4,1) ; I_phase(5,1) ; -I_phase(11,1) ; I_phase(6,1) ; -I_phase(12,1) ; -I_phase(1,1)];
It1=[It1;-It1;It1;-It1];

It2=[I_phase(1,1) ; -I_phase(8,1) ; I_phase(2,1) ; I_phase(3,1) ; -I_phase(9,1) ; I_phase(4,1) ; -I_phase(10,1) ; -I_phase(11,1) ; I_phase(5,1) ; -I_phase(12,1) ; I_phase(6,1) ; I_phase(7,1)];
It2=[It2;-It2;It2;-It2];

Jt1=It1./S_slot_down;                % 下层绕组的电流密度
Jt2=It2./S_slot_up;                  % 上层绕组的电流密度

%% permanent magnet and air-gap 永磁体与气隙之间
Ik=eye(K,K);
K_matrix=eye(K,K);
K_special=eye(K,K);
for i = 1:K
    if i == 1
        K_special(i,i)=sqrt(2);
    else
        K_matrix(i,i)=i;
        K_special(i,i)=i;
    end
end

g1=Rr/Rm;
g2=Rm/Rs;
G1=eye(K,K);
G2=eye(K,K);
for i = 1:K
    G1(i,i)=g1^i;
    G2(i,i)=g2^i;
end
K11=Ik+G1^2;
K22=Ik+G1^2;
K13=-G2;
K25=-G2;
K14=-Ik;
K26=-Ik;

K31=Ik-G1^2;
K42=Ik-G1^2;
K33=-ur*G2;
K45=-ur*G2;
K34=ur*Ik;
K46=ur*Ik;

%% 永磁体与气隙之间的Y向量

Y1=-u0*(K_special^2-Ik)^(-1)*((Rr*K*G1+Rm*Ik)*Mack-(Rr*G1+Rm*K_matrix)*Mrsk);
Y2=-u0*(K_special^2-Ik)^(-1)*((Rr*K*G1+Rm*Ik)*Mask+(Rr*G1+Rm*K_matrix)*Mrck);

Y3=-u0*(K_special^2-Ik)^(-1)*(K_matrix*(Rm*Ik-Rr*G1)*Mack-(Rm*Ik-Rr*G1)*Mrsk);
Y4=-u0*(K_special^2-Ik)^(-1)*(K_matrix*(Rm*Ik-Rr*G1)*Mask+(Rm*Ik-Rr*G1)*Mrck);

%% slot opening and slot 槽开口与槽之间
Fm=eye(M,M);
En=eye(N,N);
G3=eye(N,N);
G4=eye(M,M);
IN=eye(N,N);
gamma=eye(M,N);

for i = 1:M
    Fm(i,i)=i*pi/boa;
    G4(i,i)=(Rs/Rt)^(i*pi/boa);
end
for i = 1:N
    En(i,i)=i*pi/bsa;
    G3(i,i)=(Rt/Rsb)^(i*pi/bsa);
end
for i = 1:M
    for o = 1:N
        en=o*pi/bsa;
        fm=i*pi/boa;
        gamma(i,o)=-(2/bsa)*(en/(fm^2-en^2))*(cos(i*pi)*sin(en*(bsa+boa)/2)-sin(en*(bsa-boa)/2));
    end
end
K97=gamma'*Fm;
K98=-gamma'*Fm*G4;
K99=-En*(G3^2-IN);
K97i=eye(N*Ns,M*Ns);
K98i=eye(N*Ns,M*Ns);
K99i=eye(N*Ns,N*Ns);
Ftm=eye(M*Ns,M*Ns);%%
G4t=eye(M*Ns,M*Ns);%%
for i = 1:N
    for o = 1:M
        for n = 0:Ns-1
            K97i(i+N*n,o+M*n)=K97(i,o);
            K98i(i+N*n,o+M*n)=K98(i,o);
        end
    end
end
for i = 1:N
    for o = 1:N
        for n = 0:Ns-1
            K99i(i+N*n,o+N*n)=K99(i,o);
            Ftm(i+N*n,o+N*n)=Fm(i,o);
            G4t(i+N*n,o+N*n)=G4(i,o);
        end
    end
end

zeta=(bsa/boa)*gamma;
IM=eye(M,M);
K87i=eye(M*Ns,M*Ns);
K88i=eye(M*Ns,M*Ns);
K89i=eye(N*Ns,N*Ns);
for i = 1:M
    for o = 1:M
        for n = 0:Ns-1
            K88i(i+M*n,o+M*n)=G4(i,o);
        end
    end
end

K67=-zeta*(G3^2+IN);
for i = 1:N
    for o = 1:N
        for n = 0:Ns-1
            K89i(i+N*n,o+N*n)=K67(i,o);
        end
    end
end

Y6=zeros(N*Ns,1);
Y7=zeros(N*Ns,1);

% 双层叠绕组的Y7向量
for ns = 1:Ns
    Ji1=Jt1(ns,1);
    Ji2=Jt2(ns,1);
    for n = 1:N
        gamma0=4*cos(n*pi/2)*sin(en*boa/2)/(n*pi);
        Y7(n+(ns-1)*N,1)=-(u0/2)*(bsa/boa)*((Rsm^2-Rt^2)*Ji1+(Rsb^2-Rsm^2)*Ji2)*gamma0;
    end
end
%% slot opening and air gap 槽开口与气隙之间

K53=-K_matrix;
K65=-K_matrix;
K54=G2*K_matrix;
K66=G2*K_matrix;

yita=zeros(K,M*Ns);
kesei=zeros(K,M*Ns);
yita0=zeros(K,Ns);
kesei0=zeros(K,Ns);

for m = 1:M
    for k = 1:K
        for ns = 0:Ns-1
            fm=m*pi/boa;
            ai=theta_s/2+theta_s*ns;
            yita(k,m+M*ns)=-(1/pi)*(k/(fm^2-k^2))*(cos(m*pi)*sin(k*ai+k*boa/2)-sin(k*ai-k*boa/2));
            kesei(k,m+M*ns)=(1/pi)*(k/(fm^2-k^2))*(cos(m*pi)*cos(k*ai+k*boa/2)-cos(k*ai-k*boa/2));
        end
    end
end

for k = 1:K
    for ns = 1:Ns
        ai=theta_s/2+theta_s*(ns-1);
        yita0(k,ns)=2*sin(k*boa/2)*cos(k*ai)/(k*pi);
        kesei0(k,ns)=2*sin(k*boa/2)*sin(k*ai)/(k*pi);
    end
end

K57=yita*Ftm*G4t;
K58=-yita*Ftm;
K67=kesei*Ftm*G4t;
K68=-kesei*Ftm;

sigma=(2*pi/boa)*yita';
tao=(2*pi/boa)*kesei';
K73=sigma;
K74=sigma*G2;
K75=tao;
K76=tao*G2;
K77=-G4t;
K78=-eye(M*Ns,M*Ns);

Y8=-(u0/2)*(bsa/boa)*yita0*((Rsm^2-Rt^2)*Jt1+(Rsb^2-Rsm^2)*Jt2);
Y9=-(u0/2)*(bsa/boa)*kesei0*((Rsm^2-Rt^2)*Jt1+(Rsb^2-Rsm^2)*Jt2);


%% 拼凑系数矩阵

Z1=zeros(K,K);
Z2=zeros(K,M*Ns);
Z3=zeros(M*Ns,K);
Z4=zeros(M*Ns,M*Ns);

ratio=[K11,Z1,K13,K14,Z1,Z1,Z2,Z2,Z2;
    Z1,K22,Z1,Z1,K25,K26,Z2,Z2,Z2;
    K31,Z1,K33,K34,Z1,Z1,Z2,Z2,Z2;
    Z1,K42,Z1,Z1,K45,K46,Z2,Z2,Z2;     % 前四列的边界条件为永磁体与气隙之间
    Z1,Z1,K53,K54,Z1,Z1,K57,K58,Z2;
    Z1,Z1,Z1,Z1,K65,K66,K67,K68,Z2;
    Z3,Z3,K73,K74,K75,K76,K77,K78,Z4;
    Z3,Z3,Z3,Z3,Z3,Z3,K87i,K88i,K89i;
    Z3,Z3,Z3,Z3,Z3,Z3,K97i,K98i,K99i];

Y=[Y1;Y2;Y3;Y4;Y8;Y9;zeros(M*Ns,1);Y6;Y7];

%% 求出系数向量

R=(ratio^(-1))*Y;

A1=R(1:K,1);
C1=R(K+1:2*K,1);
A2=R(2*K+1:3*K,1);
B2=R(3*K+1:4*K,1);

C2=R(4*K+1:5*K,1);
D2=R(5*K+1:6*K,1);
C4t=R(6*K+1:6*K+M*Ns,1);
D4t=R(6*K+M*Ns+1:6*K+M*Ns*2,1);
D3t=R(6*K+M*Ns*2+1:6*K+M*Ns*3,1);


%% 求气隙磁密分布

r=0.10775;     % r是所取的气隙半径，计算的是半径为r处的径向磁密

accuracy=6000;
alpha=linspace(0,pi,accuracy);
B2r=zeros(accuracy,1);
B2a=zeros(accuracy,1);
for i = 1:accuracy
    level1=0;
    level2=0;
    level3=0;
    level4=0;
    for j = 1:ii
        n=2*j-1;
        k=p*n;
        level1=level1+k*((A2(k)/Rs)*(r/Rs)^(k-1)+(B2(k)/Rm)*(r/Rm)^(-k-1))*sin(p*k*alpha(i));
        level2=level2+k*((C2(k)/Rs)*(r/Rs)^(k-1)+(D2(k)/Rm)*(r/Rm)^(-k-1))*cos(p*k*alpha(i));
        level3=level3+k*((A2(k)/Rs)*(r/Rs)^(k-1)-(B2(k)/Rm)*(r/Rm)^(-k-1))*cos(p*k*alpha(i));
        level4=level4+k*((C2(k)/Rs)*(r/Rs)^(k-1)-(D2(k)/Rm)*(r/Rm)^(-k-1))*sin(p*k*alpha(i));
    end
    B2r(i)=-level1+level2;
    B2a(i)=-level3-level4;
end

%% 负载气隙磁密波形

filename1="D:\multi_phase_shift_angle\phase_12_Br_30_deg.csv";
filename2="D:\multi_phase_shift_angle\phase_12_Bt_30_deg.csv";
data1=readmatrix(filename1);
data2=readmatrix(filename2);
FEA_Br=data1(:,2);
FEA_Bt=data2(:,2);


figure

subplot(2,1,1)
plot(p*180*(alpha)/pi,B2r,'-.r','DisplayName',"解析法");
hold on
plot(p*180*(alpha)/pi,FEA_Br,'-k','DisplayName',"有限元法")
hold on
plot(p*180*(alpha)/pi,B2r-FEA_Br,':b','DisplayName',"绝对误差")
hold on

avg1=sum(abs(B2r-FEA_Br))/6000;
txt1=["平均绝对误差为" string(avg1)];
text(250,0.9,txt1,'FontSize',14);

xlim([0,360]);
xlabel("转子位置角度 / Degree")
ylabel("气隙磁密 / Tesla")
set(gca,'FontSize',14)
title("径向带载气隙磁密")
legend

subplot(2,1,2)
plot(p*180*(alpha)/pi,B2a,'-.r','DisplayName',"解析法");
hold on
plot(p*180*(alpha)/pi,FEA_Bt,'-k','DisplayName',"有限元法")
hold on
plot(p*180*(alpha)/pi,B2a-FEA_Bt,':b','DisplayName',"绝对误差")
hold on

avg2=sum(abs(B2a-FEA_Bt))/6000;
txt2=["平均绝对误差为" string(avg2)];
text(250,0.3,txt2,'FontSize',14);

xlim([0,360]);
xlabel("转子位置角度 / Degree")
ylabel("气隙磁密 / Tesla")
set(gca,'FontSize',14)
title("周向带载气隙磁密")
legend

%% 负载电磁力密度波形

filename3="D:\multi_phase_shift_angle\FEA_Fr_30_deg.csv";
filename4="D:\multi_phase_shift_angle\FEA_Ft_30_deg.csv";

data3=readmatrix(filename3);
data4=readmatrix(filename4);

FEA_Fr=data3(:,2);
FEA_Ft=data4(:,2);

fr=(B2r.^2-B2a.^2)/(2*u0);   % 解析法算出来的径向电磁力密度
ft=(B2r.*B2a)/u0/2;            % 解析法算出来的周向电磁力密度
% 画电磁力密度波形
figure

subplot(2,1,1)
plot(180*p*alpha/pi,fr,'-.r','DisplayName',"解析法");
hold on
plot(180*p*alpha/pi,FEA_Fr,'-k','DisplayName',"有限元法");
hold on
plot(180*p*alpha/pi,fr-FEA_Fr,':b','DisplayName',"绝对误差");
hold on
legend

avg3=sum(abs(fr-FEA_Fr))/6000;
txt3=["平均绝对误差为" string(avg3)];
text(250,5.4e5,txt3,'FontSize',14);

xlim([0 360]);
xlabel("转子位置角度 / Degree")
ylabel("径向电磁力密度 / N/m^2")
set(gca,'FontSize',14)
title("径向电磁力密度对比")

subplot(2,1,2)
plot(180*p*alpha/pi,ft,'-.r','DisplayName',"解析法");
hold on
plot(180*p*alpha/pi,FEA_Ft,'-k','DisplayName',"有限元法");
hold on
plot(180*p*alpha/pi,ft-FEA_Ft,':b','DisplayName',"绝对误差");
hold on
legend


avg4=sum(abs(ft-FEA_Ft))/6000;
txt4=["平均绝对误差为" string(avg4)];
text(250,1.65e5,txt4,'FontSize',14);

xlim([0 360]);
xlabel("转子位置角度 / Degree")
ylabel("周向电磁力密度 / N/m^2")
set(gca,'FontSize',14)
title("周向电磁力密度对比")
