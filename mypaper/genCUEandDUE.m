function [Flag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE,numDUE)
%首先，根据泊松分布（泊松分布的参数是车辆之间的平均距离）在高速公路上生成汽车。
%然后，生成DUE和CUE的用户位置

%输入：d0:高速公路的长度
%      d_avg:车辆间的平均距离，d_avg = 2.5*v
%      其他的参数可以从命名中推出
%输出：Flag，0/1，如果生成的汽车数量不够，flag=1
%     vehPos,N*2的矩阵，存储汽车坐标
%     indCUE, numCUE*1的矩阵,存储CUE用户坐标
%     indDUE, numDUEx1的矩阵，存储DUE用户发射端的坐标
%     indDUE2, numDUEx1的矩阵，存储DUE用户接收端的坐标（最接近DUE发射端的用户）
vehPos = []; % 生成汽车用户的坐标
indCUE = [];
indDUE = [];
indDUE2 = [];
Flag = 0;
%% 生成所有汽车用户的坐标并存在vehPos里面
for ilane = 1:numLane
    npoints = poissrnd(2*d0/d_avg);
    pproc = (rand(npoints,1)*2-1)*d0; % 水平坐标
    pproc2 = [pproc, (disBstoHwy+ilane*laneWidth)*ones(length(pproc), 1)]; % 地平线垂直坐标
    vehPos = [vehPos; pproc2];
end
numVeh = size(vehPos, 1);
if numVeh < numCUE+2*numDUE
    Flag = 1; % 生成的汽车用户不够多
    return; 
end
%% 生成CUE和DUE用户坐标
indPerm = randperm(numVeh);
indDUE = indPerm(1:numDUE); % randomly pick numDUE DUEs
indDUE2 = zeros(1,numDUE); % corresponding DUE receiver
for ii = 1 : numDUE
    %为每一个indDUE用户里面的元素分配一个最近的用户存在indDUE2里面
    minDist = 2*d0;
    tmpInd = 0;
    for iii = 1:numVeh
        if any(abs(iii-indDUE)<1e-6) || any(abs(iii-indDUE2)<1e-6) % iii在indDUE里面或者indDUE2里面
            continue;
        end
        newDist = sqrt((vehPos(indDUE(ii), 1)-vehPos(iii,1))^2 + (vehPos(indDUE(ii), 2)-vehPos(iii,2))^2);
        if newDist < minDist
            tmpInd = iii;
            minDist = newDist;
        end
    end
    indDUE2(ii) = tmpInd; % 最近的DUE用户对
end

cntCUE = numDUE+1; % 排除在indDUE里面的用户
while cntCUE <= numVeh
    if any(abs(indPerm(cntCUE)-indDUE2)<1e-6) % 在indDUE2里面的元素什么也不做
    else
        indCUE = [indCUE indPerm(cntCUE)];
    end
    cntCUE = cntCUE + 1;
    if length(indCUE) >= numCUE
        break
    end
end
indCUE = indCUE(:);
indDUE = indDUE(:);
indDUE2 = indDUE2(:);


%画图程序
x0=0;                          %画出圆形区域的边框
y0=0;
r=d0;
theta=0:pi/50:2*pi;
x1=x0+r*cos(theta);
y1=y0+r*sin(theta);
x=vehPos(indCUE,1);
y=vehPos(indCUE,2);
plot(x0,y0,'sk',x,y,'v');
hold on;
x2=vehPos(indDUE,1);
y2=vehPos(indDUE,2);
x3=vehPos(indDUE2,1);
y3=vehPos(indDUE2,2);
plot(x2,y2,'r*',x3,y3,'b*');
plot(x1,y1,'-k')
hold on
hold off
legend('BS','CU','DTx','DRx');
axis equal

end

