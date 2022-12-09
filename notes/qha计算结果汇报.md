## 前言
1. $VO_2$相变如何发现
2. 发现$VO_2$的两种结构
3. DFT计算得到相变温度与实验相符
4. $VO_2$材料有何特性，导致运算量上升，引出可使用+U的办法
5. +U有什么问题，U应该加多少，有什么限制
6. 本文采用什么相变计算方法，前人有谁做过，从相变能量角度应该考虑+U为多少

## 正文
### 结构
M结构(12原子）
> 本文使用jpg格式，latex中使用pdf格式

改颜色

![](D:\repeat_others_work\VO2\M\POSCAR_M.jpg "M结构")

R结构（6原子）

![](D:\repeat_others_work\VO2\R\POSCAR_R.jpg "R结构")

增加一个表格

### qha
#### 计算方法
1. 将结构都变成12原子结构
2. 对两种结构加入三组U值分别为：1，1.5，2
2. 对每一种进行`ISIF=3`的弛豫，达到稳定状态
3. 将POSCAR中的factor分别设置为：0.980， 0.985，0.990，0.995，1.005，1.010，1.015，1.020，进行`ISIF=4`的等体积弛豫，相当于在一定压强下优化
4. 使用更大的晶胞即221倍为超晶胞，计算声子谱，确定没有虚频后计算其热性质，加一个图
5. 对其进行qha计算，得到0压力下Gibbs能量随着温度的变化
6. 两种结构，三组U值，分别考虑不同的处理方法（是否使用efe），共得到12组gibbs能量，将两种结构对应相同的U值和处理方法进行比较，可以知道在某一温度相交时即对应相变发生位置

**helmholtz-volume的图**

#### 计算结果
![](D:\repeat_others_work\VO2\qha\u1.0-ev.jpg)
![](D:\repeat_others_work\VO2\qha\u2.0-ev.jpg)
![](D:\repeat_others_work\VO2\qha\u1.5-ev.jpg)
![](D:\repeat_others_work\VO2\qha\u1.0-efe.jpg)
![](D:\repeat_others_work\VO2\qha\u2.0-efe.jpg)
![](D:\repeat_others_work\VO2\qha\u1.5-efe.jpg)


标题放在里面

增加一个小图

可以看到U=1时交点过于靠前（甚至没有交点），U=2时交点过于靠后，选取U=1.5时使用efe恰好使得交点位于340K，即相变的温度。

### U的范围
以上叙述了U在考虑相变时的取值，那么U必须满足什么条件呢？
> 参考文献

我们在此基础上进行了pbe+U的计算，得到不同的能带间隔：

| U | M      | R        |
|---|--------|----------|
| 1 | 0.2138 | 0  |
| 2 | 0.4424 | 0  |
| 3 | 0.7073 | 0  |
| 4 | 0.9994 | 0  |

可以看到在增大U的时候，M不会从非金属变为金属，但是R却可能从金属变为非金属，所以U的大小有一定的限制

### 能带和声子谱
U=1.5

![](D:\repeat_others_work\VO2\M_phon_band.jpg)

![](D:\repeat_others_work\VO2\R_phon_band.jpg)

U=1.5
![](D:\repeat_others_work\VO2\M_band.jpg)

![](D:\repeat_others_work\VO2\R_band.jpg)

## 问题

1. 画的图线条以及坐标轴粗细，字体大小







