% 论文的主函数
% 比较不同车载速度下V2I链路的最大速率和最小速率

tic
clear;
clc

channNum = 1e3;
rng(3); % 控制伪随机数的产生

%% 建立参数
infty = 2000; % 用作仿真的无穷大值
dB_Pd_max = 23; % DUE用户的最大发送功率（dBm）
dB_Pc_max = 23; % CUE用户的最大发送功率（dBm）

% 大尺度衰落参数
stdV2V = 3; % 阴影衰落的标准差
stdV2I = 8;

% 蜂窝网络参数的建立
freq = 2; % 载波频率（Ghz）
radius = 500; % 蜂窝半径（m）
bsHgt = 25; % 基站高度（m）
disBstoHwy = 35; % 基站-高速公路之间的距离（m）
bsAntGain = 8; % 基站天线增益8 dBi
bsNoiseFigure = 5; % 基站噪声指数 5 dB

vehHgt = 1.5; % 汽车天线高度（m）
vehAntGain = 3; % 汽车天线增益3 dBi
vehNoiseFigure = 9; % 汽车噪声指数9 dB

numLane = 6;
laneWidth = 4;
v = 60:10:140; % 汽车速度
d_avg_ = 2.5.*v/3.6; % 根据TR366.85协议车辆间的平均距离

% V2I链路和V2V链路的QoS参数
r0 = 0.5; % CUE用户最小速率 （bps/Hz）
dB_gamma0 = 5; % DUE用户最小的SINR（dB）
p0 = 0.001; % DUE用户的中断概率阈值
dB_sig2 = -114; % 噪声功率（dBm）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% dB到线性尺度的转换
sig2 = 10^(dB_sig2/10);
gamma0 = 10^(dB_gamma0/10);
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
numCUE = 20;
numDUE = 20;
sumRate_maxSum = zeros(length(d_avg_), 1);
minRate_maxSum = zeros(length(d_avg_), 1);
sumRate_maxMin = zeros(length(d_avg_), 1);
minRate_maxMin = zeros(length(d_avg_), 1);
%%

parfor ind = 1 : length(sumRate_maxSum)
    d_avg = d_avg_(ind);
    
    cntChann = 0; % 实现信道计数
    while cntChann < channNum
        %% 在高速公路上随机产生车辆
        d0 = sqrt(radius^2-disBstoHwy^2);
        [genFlag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
        if genFlag == 1
            continue; % 产生车子的数量不够，跳入下一次循环
        end           %contiune表示结束本次循环进行下一次循环，在这里表示结束if循环，进行下一次的while语句。
                      %break表示提前结束循环，表示直接结束while循环。
        %% 生成随机大尺度衰落参数
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);
        for m = 1 : numCUE
            dist_mB = sqrt(vehPos(indCUE(m),1)^2 + vehPos(indCUE(m),2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, freq) + vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB/10);
            
            for k = 1 : numDUE
                dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE(k),1))^2 +  (vehPos(indCUE(m),2)-vehPos(indDUE(k),2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
                alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
            end
        end
        for k = 1 : numDUE
            dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 +  (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
            alpha_k_(k) = 10^(dB_alpha_k/10);
            
            dist_k = sqrt(vehPos(indDUE(k),1)^2 + vehPos(indDUE(k),2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_k, vehHgt, bsHgt, freq)+ vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB/10);
        end
        
        %% 资源分配设计-对于单个V2I-V2V用户对
        C_mk = zeros(numCUE, numCUE); % 如果numDUE < numCUE，生成虚拟DUE用户
        for m = 1 : numCUE
            alpha_mB = alpha_mB_(m);
            for k = 1 : numCUE
                
                if k > numDUE % numDUE < numCUE生成虚拟的DUE用户
                    a = Pc_max*alpha_mB/sig2;% 缺少DUE用户的情况下，CUE用户功率应该设为最大
                    C_mk(m,k) = computeCapacity(a,0); % CUE用户不受干扰的情况
                    continue;
                end
                
                alpha_k = alpha_k_(k);
                alpha_kB = alpha_kB_(k);
                alpha_mk = alpha_mk_(m,k);
                
                Pc_dmax = alpha_k*Pd_max/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(Pd_max*alpha_k))/(1-p0)-1);
                if Pc_dmax <= Pc_max
                    Pd_opt = Pd_max;
                    Pc_opt = Pc_dmax;
                else
                    %% 二分法查找Pd_cmax
                    epsi = 1e-5;
                    Pd_left = -gamma0*sig2/(alpha_k*log(1-p0)); % 线性规划图上横坐标轴的交点P_{k,min}^d
                    Pd_right = Pd_max;
                    tmpVeh = 0;
                    while Pd_right - Pd_left > epsi
                        tmpVeh = (Pd_left + Pd_right)/2;
                        if alpha_k*tmpVeh/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(tmpVeh*alpha_k))/(1-p0)-1) > Pc_max
                            Pd_right = tmpVeh;
                        else
                            Pd_left = tmpVeh;
                        end
                    end
                    Pd_cmax = tmpVeh;
                    
                    Pd_opt = Pd_cmax;
                    Pc_opt = Pc_max;
                end
                %% C_km 表示V2I链路的最大容量当第k个V2V链路和第m个V2I链路共享频谱的时候
                a = Pc_opt*alpha_mB/sig2;
                b = Pd_opt*alpha_kB/sig2;
                C_mk(m,k) = computeCapacity(a,b);
                if C_mk(m,k) < r0 % V2I链路的最小速率
                    C_mk(m,k) = -infty;
                end
            end
        end
        
        %% 复用对匹配
        [assignmentSum, cost] = munkres(-C_mk);
        [sumVal_sum, minVal_sum] = sumAndMin(C_mk, assignmentSum);
        
        [assignmentMin, dummyMin ] = maxMin( C_mk );
        [sumVal_min, minVal_min] = sumAndMin(C_mk, assignmentMin);
        
        if minVal_sum < 0 ||  minVal_min < 0 % 不可行的情况
            continue;
        end
        
        sumRate_maxSum(ind) = sumRate_maxSum(ind) + sumVal_sum;
        minRate_maxSum(ind) = minRate_maxSum(ind) + minVal_sum;
        sumRate_maxMin(ind) = sumRate_maxMin(ind) + sumVal_min;
        minRate_maxMin(ind) = minRate_maxMin(ind) + minVal_min;
    
        cntChann = cntChann + 1;
    end
    ind
    
end

sumRate_maxSum = sumRate_maxSum/channNum;
minRate_maxSum = minRate_maxSum/channNum;
sumRate_maxMin = sumRate_maxMin/channNum;
minRate_maxMin = minRate_maxMin/channNum;
%%



%% 改变不同的发射功率重新运行一遍
dB_Pd_max = 17; 
dB_Pc_max = 17; 
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);

sumRate_maxSum2 = zeros(length(d_avg_), 1);
minRate_maxSum2 = zeros(length(d_avg_), 1);
sumRate_maxMin2 = zeros(length(d_avg_), 1);
minRate_maxMin2 = zeros(length(d_avg_), 1);
%%
parfor ind = 1 : length(sumRate_maxSum)
    d_avg = d_avg_(ind);
    
    cntChann = 0; 
    while cntChann < channNum
        %% 生成汽车用户
        d0 = sqrt(radius^2-disBstoHwy^2);
        [genFlag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE,numDUE);
        if genFlag == 1
            continue; 
        end
        
        %% 生成随机大尺度衰落参数
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);
        for m = 1 : numCUE
            dist_mB = sqrt(vehPos(indCUE(m),1)^2 + vehPos(indCUE(m),2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, freq) + vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB/10);
            
            for k = 1 : numDUE
                dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE(k),1))^2 +  (vehPos(indCUE(m),2)-vehPos(indDUE(k),2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
                alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
            end
        end
        for k = 1 : numDUE
            dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 +  (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
            alpha_k_(k) = 10^(dB_alpha_k/10);
            
            dist_k = sqrt(vehPos(indDUE(k),1)^2 + vehPos(indDUE(k),2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_k, vehHgt, bsHgt, freq)+ vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB/10);
        end
        
        %% 资源分配设计-针对单个V2I-V2V用户
        C_mk = zeros(numCUE, numCUE);
        for m = 1 : numCUE
            alpha_mB = alpha_mB_(m);
            for k = 1 : numCUE
                
                if k > numDUE 
                    a = Pc_max*alpha_mB/sig2;
                    C_mk(m,k) = computeCapacity(a,0); 
                    continue;
                end
                
                alpha_k = alpha_k_(k);
                alpha_kB = alpha_kB_(k);
                alpha_mk = alpha_mk_(m,k);
                
                Pc_dmax = alpha_k*Pd_max/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(Pd_max*alpha_k))/(1-p0)-1);
                if Pc_dmax <= Pc_max
                    Pd_opt = Pd_max;
                    Pc_opt = Pc_dmax;
                else
                    %% 二分法找Pd_cmax
                    epsi = 1e-5;
                    Pd_left = -gamma0*sig2/(alpha_k*log(1-p0)); 
                    Pd_right = Pd_max;
                    tmpVeh = 0;
                    while Pd_right - Pd_left > epsi
                        tmpVeh = (Pd_left + Pd_right)/2;
                        if alpha_k*tmpVeh/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(tmpVeh*alpha_k))/(1-p0)-1) > Pc_max
                            Pd_right = tmpVeh;
                        else
                            Pd_left = tmpVeh;
                        end
                    end
                    Pd_cmax = tmpVeh;
                    
                    Pd_opt = Pd_cmax;
                    Pc_opt = Pc_max;
                end
                %% C_km表示当第m个V2I链路复用第k个V2V链路频谱时V2I链路的总速率
                a = Pc_opt*alpha_mB/sig2;
                b = Pd_opt*alpha_kB/sig2;
                C_mk(m,k) = computeCapacity(a,b);
                if C_mk(m,k) < r0 
                    C_mk(m,k) = -infty;
                end
            end
        end
        
        %% 复用对匹配
        [assignmentSum, cost] = munkres(-C_mk);
        [sumVal_sum, minVal_sum] = sumAndMin(C_mk, assignmentSum);
        
        [assignmentMin, dummyMin ] = maxMin( C_mk );
        [sumVal_min, minVal_min] = sumAndMin(C_mk, assignmentMin);
        
        if minVal_sum < 0 ||  minVal_min < 0 
            continue;
        end
        
        sumRate_maxSum2(ind) = sumRate_maxSum2(ind) + sumVal_sum;
        minRate_maxSum2(ind) = minRate_maxSum2(ind) + minVal_sum;
        sumRate_maxMin2(ind) = sumRate_maxMin2(ind) + sumVal_min;
        minRate_maxMin2(ind) = minRate_maxMin2(ind) + minVal_min;
    
        cntChann = cntChann + 1;
    end
    ind
end

sumRate_maxSum2 = sumRate_maxSum2/channNum;
minRate_maxSum2 = minRate_maxSum2/channNum;
sumRate_maxMin2 = sumRate_maxMin2/channNum;
minRate_maxMin2 = minRate_maxMin2/channNum;


%%画图
LineWidth = 1.5;
MarkerSize = 6;
figure
plot(v, sumRate_maxSum, 'k-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxMin, 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxSum2, 'k--s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxMin2, 'b--o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
grid on
legend('P_{max}^c = 23 dBm, Algorithm 1', 'P_{max}^c = 23 dBm, Algorithm 2',...
'P_{max}^c = 17 dBm, Algorithm 1', 'P_{max}^c = 17 dBm, Algorithm 2')
xlabel('$v$ (km/h)', 'interpreter','latex')
ylabel('$\sum\limits_m C_m$ (bps/Hz)', 'interpreter','latex')
% saveas(gcf, sprintf('sumRateVsSpeed')); % save current figure to file

figure
plot(v, minRate_maxSum, 'k-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxMin, 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxSum2, 'k--s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxMin2, 'b--o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
grid on
legend('P_{max}^c = 23 dBm, Algorithm 1', 'P_{max}^c = 23 dBm, Algorithm 2',...
'P_{max}^c = 17 dBm, Algorithm 1', 'P_{max}^c = 17 dBm, Algorithm 2')
xlabel('$v$ (km/h)', 'interpreter','latex')
ylabel('$\min C_m$ (bps/Hz)', 'interpreter','latex')
% saveas(gcf, 'minRateVsSpeed'); % save current figure to file

% save all data
save('main_rateSpeed')

toc