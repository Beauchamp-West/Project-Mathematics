#! /bin/bash
#
gfortran -c -Wall task1.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran task1.o mgmres.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm task1.o
#
mv a.out task1
# ./task1 > task1.txt
./task1
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm task1
#
echo "Normal end of execution."