#include "header.h"

PS::matrix CalcTensorProduct(const PS::F64vec arg1, const PS::F64vec arg2){
  PS::matrix tensor;
  tensor.xx = arg1.x * arg2.x;
  tensor.xy = arg1.x * arg2.y;
  tensor.xz = arg1.x * arg2.z;
  tensor.yx = arg1.y * arg2.x;
  tensor.yy = arg1.y * arg2.y;
  tensor.yz = arg1.y * arg2.z;
  tensor.zx = arg1.z * arg2.x;
  tensor.zy = arg1.z * arg2.y;
  tensor.zz = arg1.z * arg2.z;
  
  return(tensor);
}

PS::F64 CalcTensorToScalar(const PS::matrix arg1, const PS::matrix arg2){
  PS::F64 scalar;
  scalar = arg1.xx * arg2.xx + arg1.xy * arg2.xy + arg1.xz * arg2.xz + arg1.yx * arg2.yx + arg1.yy * arg2.yy + arg1.yz * arg2.yz + arg1.zx * arg2.zx + arg1.zy * arg2.zy + arg1.zz * arg2.zz;
  return(scalar);
}
