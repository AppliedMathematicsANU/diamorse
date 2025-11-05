#!/bin/bash

set -o errexit

bindir=${BINDIR:-$(dirname $0)/../bin}
timeCmd=${TIME:-"/usr/bin/time -v"}

input=$(echo $1 | sed 's/[\._]nc$//')
name=$(basename $input)
shift 1

rm -f $log

if [[ $name =~ ^segmented ]]
then
  echo "=== SignedEuclideanDistanceTransform ===" >&2
  scalars=$(echo $name | sed 's/^segmented/tomo_float/')_SEDT
  $timeCmd $bindir/SEDT ${input}[._]nc ${scalars}.nc
  echo >&2
else
  scalars=$input
fi

field=$(basename $scalars | sed 's/^tomo_float/vector_field/')_GVF_SMP

echo "=== VectorField ===" >&2
$timeCmd $bindir/VectorField ${scalars}[._]nc __field.nc
echo >&2

echo '=== Simplify ===' >&2
$timeCmd $bindir/Simplify ${scalars}[._]nc __field[._]nc ${field}.nc $*
echo >&2

pairs=$(echo $field | sed 's/^vector_field/persistence/')_PP

echo "=== Persistence ===" >&2
$timeCmd $bindir/PersistencePairs ${scalars}[._]nc ${field}[._]nc ${pairs}.txt
echo >&2

rm -rf __*
