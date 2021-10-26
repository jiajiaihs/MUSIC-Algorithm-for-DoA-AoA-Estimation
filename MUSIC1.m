clear all
close all
derad = pi/180; %角度->弧度
radeg = 180/pi; %弧度->角度
twpi = 2*pi;
kelm = 8;       %阵元个数     
dd = 0.5;       %阵元间距       
d=0:dd:(kelm-1)*dd;     
iwave = 3;      %信源数             
theta = [15 28 60];  %波达方向/待估计角度为15 28 60
%这个程序里面待估计角度是设置好的。实际应用中不应该是不知道角度去估计这个角度么？那如何生成导向矩阵A？
%实际不用，接收的数据直接为X(t)，公式中是将接收数据分解为A*S+N，所以才需要生成导向矩阵。
%导向矢量不用去生成。这是接收天线本身的属性，实际中我们只要输出。仿真中需要生成导向矢量是因为我们需要模拟天线接收到的信号，所以会给出角度。
snr = 10;            %信噪比
n = 500;             %快拍数或者称为采样数
A=exp(-1i*twpi*d.'*sin(theta*derad));    %方向向量,.'转置，'共轭转置
S=randn(iwave,n);    %信源信号，正态分布随机矩阵
X=A*S;               %接收信号
X1=awgn(X,snr,'measured');  %添加高斯白噪声
Rxx=X1*X1'/n;        %计算协方差矩阵
[EV,D]=eig(Rxx);     %计算Rxx的特征值对应的对角阵D和特征向量构成的矩阵EV
EVA=diag(D)';        %diag抽取矩阵对角线元素
[EVA,I]=sort(EVA);   %特征值从小到大排序，sort只能从小到大排列
EVA=fliplr(EVA);     %特征值左右翻转，从大到小序
EV=fliplr(EV(:,I));  %对应特征向量排序

%构造MUSIC谱函数
for iang = 1:361
        angle(iang)=(iang-181)/2;
        phim=derad*angle(iang);
        a=exp(-1i*twpi*d*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);     %得到噪声子空间
        SP(iang)=(a'*a)/(a'*En*En'*a);
end

%作图
SP=abs(SP);
SPmax=max(SP);
SP=10*log10(SP/SPmax);
h=plot(angle,SP);
set(h,'Linewidth',2)
xlabel('angle (degree)')
ylabel('magnitude (dB)')
axis([-90 90 -60 0])
set(gca, 'XTick',(-90:10:90))
grid on  
