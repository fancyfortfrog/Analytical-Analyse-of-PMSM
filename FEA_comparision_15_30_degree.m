deg_15_Br="D:\multi_phase_shift_angle\phase_12_Br_15_deg.csv";
deg_15_Bt="D:\multi_phase_shift_angle\phase_12_Bt_15_deg.csv";
deg_15_Fr="D:\multi_phase_shift_angle\FEA_Fr_15_deg.csv";
deg_15_Ft="D:\multi_phase_shift_angle\FEA_Ft_15_deg.csv";

deg_30_Br="D:\multi_phase_shift_angle\phase_12_Br_30_deg.csv";
deg_30_Bt="D:\multi_phase_shift_angle\phase_12_Bt_30_deg.csv";
deg_30_Fr="D:\multi_phase_shift_angle\FEA_Fr_30_deg.csv";
deg_30_Ft="D:\multi_phase_shift_angle\FEA_Ft_30_deg.csv";

FEA_15_Br=readmatrix(deg_15_Br);
FEA_15_Br=FEA_15_Br(:,3);

FEA_15_Bt=readmatrix(deg_15_Bt);
FEA_15_Bt=FEA_15_Bt(:,3);

FEA_15_Fr=readmatrix(deg_15_Fr);
FEA_15_Fr=FEA_15_Fr(:,3);

FEA_15_Ft=readmatrix(deg_15_Ft);
FEA_15_Ft=FEA_15_Ft(:,3);

FEA_30_Br=readmatrix(deg_30_Br);
FEA_30_Br=FEA_30_Br(:,3);

FEA_30_Bt=readmatrix(deg_30_Bt);
FEA_30_Bt=FEA_30_Bt(:,3);

FEA_30_Fr=readmatrix(deg_30_Fr);
FEA_30_Fr=FEA_30_Fr(:,3);

FEA_30_Ft=readmatrix(deg_30_Ft);
FEA_30_Ft=FEA_30_Ft(:,3);

accuracy=6000;
alpha=linspace(0,pi,accuracy);
p=2;

figure
plot(180*p*alpha/pi,FEA_15_Br,'-.r','DisplayName',"PSA=15 deg");
hold on
plot(180*p*alpha/pi,FEA_30_Br,'-k','DisplayName',"PSA=30 deg");
hold on

title("径向气隙磁密对比")
xlabel("空间角度 / DEG")
ylabel("气隙磁密 / T")
xlim([0 360])
set(gca,'FontSize',14)
legend

figure
plot(180*p*alpha/pi,FEA_15_Bt,'-.r','DisplayName',"PSA=15 deg");
hold on
plot(180*p*alpha/pi,FEA_30_Bt,'-k','DisplayName',"PSA=30 deg");
hold on

title("周向气隙磁密对比")
xlabel("空间角度 / DEG")
ylabel("气隙磁密 / T")
xlim([0 360])
set(gca,'FontSize',14)
legend

figure
plot(180*p*alpha/pi,FEA_15_Fr,'-.r','DisplayName',"PSA=15 deg");
hold on
plot(180*p*alpha/pi,FEA_30_Fr,'-k','DisplayName',"PSA=30 deg");
hold on

title("径向电磁力密度对比")
xlabel("空间角度 / DEG")
ylabel("电磁力密度 / N/m^2")
xlim([0 360])
set(gca,'FontSize',14)
legend

figure
plot(180*p*alpha/pi,FEA_15_Ft,'-.r','DisplayName',"PSA=15 deg");
hold on
plot(180*p*alpha/pi,FEA_30_Ft,'-k','DisplayName',"PSA=30 deg");
hold on

title("周向电磁力密度对比")
xlabel("空间角度 / DEG")
ylabel("电磁力密度 / N/m^2")
xlim([0 360])
set(gca,'FontSize',14)
legend