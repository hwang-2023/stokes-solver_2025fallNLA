函数用途
data.m: 用于对不同参数进行求解, 输出图表
solve_stokes.m: 对于一组给定的参数进行求解
initialize.m: 初始化$F,G$及计算真解 
left_multiply_matrix.m: 对矩阵$U,V$的向量形式左乘$[A,B;B^T,0]$后化回矩阵形式
left_multiply_A.m: 对矩阵$U,V$的向量形式左乘$A$后化回矩阵形式
left_multiply_B.m: 对矩阵$P$的向量形式左乘$B$后化回矩阵形式
left_multiply_BT.m: 对矩阵$U,V$的向量形式左乘$B^T$后化回矩阵形式
dot_multiply.m: 对两个矩阵, 化成向量形式后求内积
Uzawa.m: 实现Uzawa迭代
V_cycle.m: 实现对方程的 V-cycle方法求解 


参数设置
最大Uzawa迭代次数: MAXITER=20


如何复现上机报告中的结果
对于一组给定参数进行求解: 在solve_stokes.m中设置对应参数并运行
输出上机报告中图表: 直接运行data.m