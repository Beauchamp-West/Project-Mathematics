#! /bin/bash
#
gfortran -c -Wall pde1.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran pde1.o mgmres.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pde1.o
#
mv a.out pde1
./pde1 > pde1.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm pde1
#
echo "Normal end of execution."