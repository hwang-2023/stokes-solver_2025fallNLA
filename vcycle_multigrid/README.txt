函数用途
data.m: 用于对不同参数进行求解, 输出图表
solve_stokes.m: 对于一组给定的参数进行求解
initialize.m: 初始化$F,G$及计算真解 
left_multiply_matrix.m: 对矩阵$U,V$的向量形式左乘$[A,B;B^T,0]$后化回矩阵形式
limit.m: 对$F,G,Z$应用限制算子, 保持其矩阵形式
ift.m: 对$U,V,P$应用提升算子, 保持其矩阵形式 
DGS.m: 实现GDS迭代
V_cycle.m: 实现对方程的 V-cycle方法求解 


参数设置
V-cycle中上升, 下降分别的迭代次数: v1=v2=2 or 4 or 8
V-cycle中底层的大小: L=2 or 4
最大V-cycle次数: MAXITER=20


如何复现上机报告中的结果
对于一组给定参数进行求解: 在solve_stokes.m中设置对应参数并运行
输出上机报告中图表: 直接运行data.m