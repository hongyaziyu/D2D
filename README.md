# vehicular communication的文章部分源代码

  1. main_rateVsSpeed是主函数文件，可以显示出汽车速度对汽车通信遍历容量的影响。
  2.首先根据genCUEandDUE.m文件根据泊松分布随机生成汽车位置。然后根据genPL.m文件生成大尺度信道信道衰落参数。最后为了解决模型提出的问题，
我们首先用computeCapacity.m文件利用二分法解决非线性规划的问题。然后用munkres.m文件中匈牙利算法解决了最大权值的二分匹配问题。之后发现匈牙利算法
只能实现遍历容量最大，并没有考虑遍历容量分配的公平性问题，基于公平性角度出发，对匈牙利算法进行了一些改进，在maxMin.m文件中可以看出匈牙利算法的改进。
