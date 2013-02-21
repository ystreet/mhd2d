# compile and link script for mhd2d code
#
base_pth="/home/mdw/mhd2d"
exec_name="mhd2d_run"
#
obj_dir="$base_pth"
#echo " obj_dir = $obj_dir"
mod_dir="$base_pth"
#echo " mod_dir = $mod_dir"
src_dir="$base_pth"
#echo " src dir = $src_dir"
debug_opts="-g -traceback -check bounds"

# Lapack
f95_lapack_dir="/opt/intel/Compiler/11.1/073/mkl/lib/em64t"
f95_lapack_mod="/opt/intel/Compiler/11.1/073/mkl/include/em64t/lp64"
lnk_lapack="-zero -save -lmkl_lapack95_lp64 -mkl=sequential"
# blas_dir="/usr/lib64"
# lnk_lapack="$lapack_dir/liblapack.a -L$blas_dir -l:libblas.so.3"

# COMPILE commands and debug switches
#
echo "compiling mhd2d_constants.f90..."
ifort -c ${debug_opts} ${src_dir}/mhd2d_constants.f90
#mv mhd2d_constants.o ${obj_dir}
#mv mhd2d_constants.mod ${mod_dir}

echo "compiling mhd2d_va_cond.f90..."
ifort -c ${debug_opts} ${src_dir}/mhd2d_va_cond.f90
#mv mhd2d_va_cond.o ${obj_dir}

echo "compiling mhd2d_grid.f90..."
ifort -c ${debug_opts} ${src_dir}/mhd2d_grid.f90
#mv mhd2d_grid.o ${obj_dir}
#
echo "compiling mhd2d.f90..."
ifort -c ${debug_opts} ${src_dir}/mhd2d.f90
#mv mhd2d.o ${obj_dir}
#
# Link all these with libraries
#
echo "Linking..."
#gfortran ${f90flags} ${obj_dir}/mhd2d_va_cond.o ${obj_dir}/mhd2d_constants.o ${obj_dir}/mhd2d_grid.o ${obj_dir}/mhd2d.o -o ${exec_name} ${lnk_lapack} ${lnk_aacgm} ${debug_opts}
ifort ${f90flags} ${obj_dir}/mhd2d_va_cond.o ${obj_dir}/mhd2d_constants.o ${obj_dir}/mhd2d_grid.o ${obj_dir}/mhd2d.o -o ${exec_name} ${lnk_lapack} ${debug_opts}
