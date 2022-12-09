[phonopy官网](http://phonopy.github.io/phonopy/vasp-dfpt.html#vasp-dfpt-interface)

（但是照着官网做不出来也是真的）

## 结构驰豫
声子谱计算也需要结构驰豫

确定最小能量对应晶胞后，先在phonopy中使用命令：
```
phonopy --symmetry -c POSCAR
```

得到PPOSCAR，这个是原胞，然后使用扩胞的办法（官网上说需要扩展到10埃以上）：
```
phonopy -d --dim="2 2 2" -c PPOSCAR
```
将产生的SPOSCAR这个作为之后静态计算的POSCAR。

## DFPT计算
（另一种办法是有限位移，和这个类似，可以参见phono3py计算热导率，操作基本一致）

### KPOINTS
```
Automatic mesh
0              ! number of k-points = 0 -> automatic generation scheme
Gamma          ! generate a Gamma centered grid
7 7 7          ! subdivisions N_1, N_2 and N_3 along recipr. latt. vectors
0. 0. 0.       ! optional shift of the mesh (s_1, s_2, s_3)
```
### INCAR
之后所有的#注释部分为考虑born有效电荷近似。
```
PREC = Accurate
 IBRION = 8
  EDIFF = 1.0e-08
  IALGO = 38
 ISMEAR = 0; SIGMA = 0.05
  LREAL = .FALSE.
ADDGRID = .TRUE.
  LWAVE = .FALSE.
 LCHARG = .FALSE.
LEPSILON = .TRUE.
```
## 产生FORCE_CONSTANTS

首先处理vasprun.xml得到FORCE_CONSTANTS：
```
phonopy --fc vasprun.xml
```
## 产生born有效电荷文件
```
phonopy-vasp-born
```
将输出内容写在born文件里面。

## 计算声子谱
自己写band.conf：
```
ATOM_NAME = Cd Te
DIM = 2 2 2
CELL_FILENAME = PPOSCAR
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
# NAC = .TRUE.
BAND = 0.0 0.0 0.0 0.5 0.0 0.5 0.375 0.375 0.75 0.0 0.0 0.0 0.5 0.5 0.5
BAND_POINTS=101
BAND_LABELS =
BAND_CONNECTION = .TRUE.
FORCE_CONSTANTS = READ
```
然后使用命令生成声子谱图：
```
phonopy band.conf -p -s
```
基本就行了

可以使用如下命令将声子谱的数据输出出来：
```
phonopy-bandplot --gnuplot > phonon.out
```
也可以通过读取××.yaml得到相应数据。

## 计算热性质

写入mesh.conf文件：
```
ATOM_NAME = Cd Te
CELL_FILENAME = PPOSCAR
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
DIM = 2 2 2
MP = 8 8 8
FORCE_CONSTANTS = READ
```

然后使用如下命令绘图：
```
phonopy -t mesh.conf -p -s
```
官网上似乎说明，```thermal_properties.yaml```中的内能和自由能应该仅仅指声子的内能和自由能。