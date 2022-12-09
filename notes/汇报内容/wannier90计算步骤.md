[toc]
## PBE计算步骤
### 静态计算
**POSCAR**
```bash
CdTe
3.3191265969887582
     0.0000000000000000    1.0000000000000000    1.0000000000000000
     1.0000000000000000    0.0000000000000000    1.0000000000000000
     1.0000000000000000    1.0000000000000000    0.0000000000000000
Cd Te
   1    1
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.2499999999999998  0.2499999999999999  0.2499999999999999
```
**KPOINTS**
```bash
Automatic mesh
0              ! number of k-points = 0 -> automatic generation scheme
Gamma          ! generate a Gamma centered grid
7 7 7          ! subdivisions N_1, N_2 and N_3 along recipr. latt. vectors
0. 0. 0.       ! optional shift of the mesh (s_1, s_2, s_3)
```
**INCAR**
```bash
PREC = Normal    ! standard precision

IBRION = -1

EDIFF = 1.0e-06  ! pretty small

ISMEAR = 0       ! always well
SIGMA = 0.05     ! maybe well
```
**job**
```bash
#! /bin/bash
#PBS -m e
#PBS -l nodes=1:ppn=12
#PBS-l walltime=12:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodes
BIN=~/544-bin/vasp_std-wann90v2
mpirun -np 12 -machinefile nodes $BIN
```
**wannier90.win**
```bash
num_bands =    24  ! set to NBANDS by VASP
num_wann = 10

dis_num_iter=200
num_iter=100
iprint=3

guiding_centres=true

begin projections
Cd : s; d
Te : s; p
end projections

bands_plot =.true.
fermi_energy = 1.9373

begin kpoint_path
L 0.5 0.5 0.5 G 0.0 0.0 0.0
G 0.0 0.0 0.0 X 0.5 0.0 0.5
end kpoint_path
```
### wannier生成数据
```bash
wannier90v2.x wannier90
```
### 绘图
```bash
gnuplot
load "wannier90_band.gnu"
```
![](D://repeat_others_work//wannier90//wannier_PBE.jpg "wannier_PBE")
### 原数据图
![](D://repeat_others_work//CdTe//PBE//band//band.jpg "PBE_band")
