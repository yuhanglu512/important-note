### 结构弛豫
### 计算0K不同体积对应能量$ F_{elec}(V,0) $
> 此处不可扩胞

#### 计算job
```bash
#!/bin/bash
#PBS -m e
#PBS -l nodes=1:ppn=12
#PBS -l walltime=54:00:00
cd $PBS_O_WORKDIR


for i in {-5..5}
do
mkdir POSCAR-$i
cp POTCAR INCAR KPOINTS POSCAR-$i/
cd POSCAR-$i


j=`echo "scale=3;1+$i/100"|bc`
cat > POSCAR <<!
cif2pos.py
   $j
     5.4052174219053342   -0.0006264332341809   -0.0004300233278269
     3.0320747641875418    4.4746953282793926   -0.0004300233236884
     3.0320747641875676    1.6078276194948409    4.1758578532430937
   Ga   O
     4     6
Direct
  0.8558811160186862  0.8558811160186868  0.8558811160186868
  0.6441188839813138  0.6441188839813132  0.6441188839813132
  0.1441188839813138  0.1441188839813132  0.1441188839813132
  0.3558811160186862  0.3558811160186868  0.3558811160186868
  0.7500000000000000  0.4469494189975330  0.0530505810024670
  0.4469494189975335  0.0530505810024669  0.7499999999999999
  0.0530505810024665  0.7500000000000001  0.4469494189975331
  0.2500000000000000  0.5530505810024670  0.9469494189975330
  0.9469494189975335  0.2499999999999999  0.5530505810024668
  0.5530505810024665  0.9469494189975332  0.2500000000000001
!

cat $PBS_NODEFILE > nodes
BIN=~/544-bin/vasp_std
mpirun -np 12 -machinefile nodes $BIN
cd ..
done
```
#### 复制到本地copy.sh
```bash
#! /bin/bash
password=
path=Ga2O3/alpha/phonopy
for i in {-5..5}
do
mkdir POSCAR-$i
cd POSCAR-$i
sshpass -p $password scp yhlu@192.168.129.71:~/$path/POSCAR-$i/POSCAR ./
sshpass -p $password scp yhlu@192.168.129.71:~/$path/POSCAR-$i/vasprun.xml ./
cd ..
done
```
### 批处理
#### mesh.conf
```bash
ATOM_NAME = Ga O
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
DIM = 1 1 1
MP = 8 8 8
NAC = .TRUE.
FORCE_CONSTANTS = READ
```
#### process.sh
```bash
#! /bin/bash


for i in {-5..5}
do
cd POSCAR-$i
phonopy --fc vasprun.xml > output
phonopy-vasp-born > born
phonopy -t ../mesh.conf >> output
cp thermal_properties.yaml ../thermal_properties.yaml-$i
cd ..
done
phonopy-vasp-efe POSCAR-{-5..5}/vasprun.xml > output
phonopy-qha -p -s --efe fe-v.dat e-v.dat thermal_properties.yaml-{-5..5} > parameters.dat
grep Parameter helmholtz-volume.dat > paras.dat
```
其中`phonopy-vasp-efe POSCAR-{-5..5}/vasprun.xml `将产生e-v.dat和fe-v.dat。
> e-v.dat是不同体积在0K下对应的内能$ F_{elec}(V,0K) $（或者说自由能），通过在OUTCAR中查询得到。
> fe-v.dat是不同体积不同温度对应的电子自由能$ F_{elec}(V,T) $，与0K时的自由能变化非常小。


其中`phonopy-qha -p -s --efe fe-v.dat e-v.dat thermal_properties.yaml-{-5..5} `将产生最终需要的所有数据。
> 1. helmholtz-volume.dat中数据为$ F_{all}(V,T) $是在fe-v.dat基础上加上声子能量$ F_{vib}(V,T) $得到的，其中声子能量由thermal_properties.yaml得到，声子能量须化为eV/unit为单位。通过计算可以证明。通过计算表明高温时的$ F_{all}(V,T)=F_{elec}(V,T)+F_{vib}(V,T) $。
> 2. parameters.dat和paras.dat中的拟合参数有3组都是相同的，但是体积模量是不同的，差了160倍，为能量换算问题，即$ eV=1.6 \times 10^{-19} \times 10^{-9} GPa \times 10^{30} \overset{\circ}{A}^3 = 160 GPa \, \overset{\circ}{A}^3 $。
> 3. helmholtz-volume_fitted.dat对总自由能做了处理，为了画图方便，但是不能直接使用。

注：
> 1. 虽然相图是关于p,T的，但是p是由V决定的
