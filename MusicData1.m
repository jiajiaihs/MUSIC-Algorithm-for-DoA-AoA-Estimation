clear all
close all
clc
derad = pi/180; % 角度->弧度
radeg = 180/pi; % 弧度->角度
twpi = 2*pi;
kelm = 8;       % 阵元个数   
dd = 0.5;       % 阵元间距   
d=0:dd:(kelm-1)*dd;
iwave = 3;      % 信源数
theta = -60:0.5:60;  % MUSIC谱峰搜索的范围
snr = 10;            % 信噪比
n = 512;             % 快拍数或者称为采样数
A=exp(-1i*twpi*d.'*sin(theta*derad)); % 方向矢量,方向向量,.'转置，'共轭转置
iq = dlmread('data1.txt',',',1,0); 
S = 1i*iq(1:2:end) + iq(2:2:end);     % 这样就有虚部实部了
X = reshape(S,kelm,n); 
% X1 = awgn(X,snr,'measured');% 在添加噪声前计算信号X的功率(dBW)
X1 = awgn(X,snr,'measured');
Rxx = X1*X1'/n;
InvS = inv(Rxx);
[EV,D] = eig(Rxx);            % 特征值分解
EVA = diag(D)';               % 将特征值矩阵对角线提取并转为一行
[EVA,I] = sort(EVA);          % 将特征值排序，从小到大
EVA = fliplr(EVA);
EV = fliplr(EV(:,I));         % 对应特征矢量排序

% 遍历每个角度，计算空间谱
for iang = 1:361
    angle(iang) = (iang-181)/2;    % 相位角
    phim = derad*angle(iang);
    a = exp(-1i*twpi*d*sin(phim)).';  
    L = iwave;
    En=EV(:,L+1:kelm);             % 得到噪声子空间
    SP(iang) = (a'*a)/(a'*En*En'*a);
end 

% 作图
SP = abs(SP);
SPmax = max(SP);             % 求最大值函数
SP = 10*log10(SP/SPmax);     % 归一化处理
h = plot(angle,SP);
set(h,'Linewidth',2)
xlabel('入射角angle (degree)')
ylabel('空间谱magnitude (dB)')
axis([-90 90 -60 0])
set(gca, 'XTick',(-90:10:90))
grid on     