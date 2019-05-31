function [ combinedPL ] = genPL(linkType, stdShadow, dist, hgtTX, hgtRX, freq)
%genPL:计算大尺度衰落参数（路径衰落+阴影效应）
%输入：linkType：选择模式，‘V2V’或者‘V2I’
%      stdShadow：对数正态分布的阴影衰落平方差（dB）
%      dist：TX和RX之间的距离（m）
%      heightTX：天线发送端的高度（m）
%      heightRX：天线接收端的高度（m）
%      freq：载波频率（Ghz）
%
%输出：combinedPL，联合大尺度衰落值（路径衰落+阴影效应）


if strcmp(upper(linkType), 'V2V')
    d_bp = 4*(hgtTX-1)*(hgtRX-1)*freq*10^9/(3*10^8);
    A = 22.7; B = 41.0; C = 20;
    if dist <= 3
        PL = A*log10(3) + B + C*log10(freq/5);
    elseif dist <= d_bp
        PL = A*log10(dist) + B + C*log10(freq/5);
    else
        PL = 40*log10(dist)+9.45-17.3*log10((hgtTX-1)*(hgtRX-1))+2.7*log10(freq/5);
    end
else
    PL = 128.1 + 37.6*log10(sqrt((hgtTX-hgtRX)^2+dist^2)/1000);
end

combinedPL = - (randn(1)*stdShadow + PL);

end