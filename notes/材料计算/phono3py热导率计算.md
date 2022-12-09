## 构建超晶胞
### 二阶与三阶相同

注意：扩胞之前最好是晶胞，似乎需要扩大到64个原子可以计算出正确的结果。

```
phono3py -d --dim="2 2 2" -c PPOSCAR
```
晶胞大就不扩胞了，小的话就扩一下胞

### 二阶与三阶不同
```
phono3py -d --dim-fc2="3 3 3" --dim="2 2 2" -c PPOSCAR
```
这个会生成POSCAR_FC2-×××××

## 将本地文件复制到集群
```
scp POSCAR* yhlu@192.168.129.71:~/phono3py/
```
## 投入vasp运算
### job
```bash
#! /bin/bash
#PBS -m e
#PBS -l nodes=1:ppn=28
#PBS -l walltime=96:00:00
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

source /etc/profile
BIN=/opt/software/vasp/5.4.4/vasp_std

for j in {00001..00114}
do
i=`printf "%05d" $j`
mkdir $i
cp POSCAR-$i $i/POSCAR
cp KPOINTS POTCAR INCAR $i/
cd $i

cat $PBS_NODEFILE > nodes
mpirun -np 28 -machinefile nodes $BIN
cd ..
done
```
如果三阶和二阶不一样，还要在后面加一段：
```
for j in {00001..00002}
do
i=`printf "%05d" $j`
mkdir FC2-$i
cp POSCAR_FC2-$i FC2-$i/POSCAR
cp POTCAR INCAR FC2-$i/
cd FC2-$i
cp ../KPOINTS_FC ./KPOINTS

cat $PBS_NODEFILE > nodes
BIN=/opt/software/vasp/5.4.4/vasp_std
mpirun -np 28 -machinefile nodes $BIN

cd ..
done
```

### KPOINTS
```
Automatic mesh
0              ! number of k-points = 0 -> automatic generation scheme
Gamma          ! generate a Gamma centered grid
7 7 7          ! subdivisions N_1, N_2 and N_3 along recipr. latt. vectors
0. 0. 0.       ! optional shift of the mesh (s_1, s_2, s_3)
```
### INCAR
不知道对不对

目前得到的结果是ADDGRID有时不影响计算速度，IALGO似乎是默认值，似乎EDIFF=1E-06也可以，并且ADDGRID不用加，至少对于Si是成立的。
```
PREC = Accurate
IBRION = -1
EDIFF = 1.0e-08
ISMEAR = 0; SIGMA = 0.05
IALGO = 38
LWAVE = .FALSE.
LCHARG = .FALSE.
ADDGRID = .TRUE.
```
## 将大量vasprun.xml复制到本地
将以下内容写在copy.sh中，然后运行bash copy.sh
```bash
#! /bin/bash

your_password=

for j in {1..114}
do
i=`printf "%05d" $j`
mkdir $i
sshpass -p $your_password scp yhlu@192.168.129.71:~/phono3py/$i/vasprun.xml $i/
done

for j in {00001..00002}
do
i=`printf "%05d" $j`
mkdir FC2-$i
sshpass -p $your_password scp yhlu@192.168.129.71:~/phono3py/FC2-$i/vasprun.xml FC2-$i/
done
```

如果三阶和二阶一样的话，是没有后一段代码的，然后就会报错，不过无所谓，之后的内容正常运行
## 产生FORCE_FC
### 二阶和三阶相同
```
phono3py --cf3 {00001..00114}/vasprun.xml
```

### 二阶和三阶不同
增加以下命令
```
phono3py --cf2 FC2-{00001..00002}/vasprun.xml
```
## 产生fc.hdf5

```bash
phono3py --sym-fc
```
计算热导率

```bash
phono3py --mesh="11 11 11" --fc3 --fc2 --br --nac -o nac
```

在此基础上还可以加入`--isotope`同位素影响，和`--bmfp=1000`（单位为微米）的影响。