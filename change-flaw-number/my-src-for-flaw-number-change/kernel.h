#pragma once

struct kernel_t{
  kernel_t(){
  }
  //W
  PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
    PS::F64 value;

    if(PARAM::KERNEL_FUNCTION==0){
      //gaussian kernel
      value = pow( h * sqrt(M_PI) , - PARAM::Dim ) * exp( - ( dr * dr ) / ( h * h ) );
    }
    if(PARAM::KERNEL_FUNCTION==1){
      //cubic spline kernel
      PS::F64 alpha = 1.0/(M_PI*h*h*h);
      PS::F64 q = sqrt(dr*dr)/h;
      if(q < 1.0){
        value = alpha*( 1.0 - 1.5*q*q + 0.75*q*q*q );
      }
      else if(q < 2.0){
        value = alpha*0.25*math::pow3(2.0-q);
      }
      else{
        value = 0.0;
      }
    }
    return value;
  }
  //gradW
  PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
    PS::F64vec vec_value;

    if(PARAM::KERNEL_FUNCTION==0){
      //gradient of gaussian kernel
      vec_value = dr * ( pow( h * sqrt(M_PI) , - PARAM::Dim ) * ( -2.0 / ( h * h ) ) * exp( - ( dr * dr ) / ( h * h ) ) );
    }
    if(PARAM::KERNEL_FUNCTION==1){
      //gradient of cubic spline kernel
      PS::F64 alpha = 9.0/(4.0*M_PI*pow(h,5.0));
      PS::F64 q = sqrt(dr*dr)/h;
      if(q < 1.0){
        vec_value = alpha*(q-1.333333333333)*dr;
      }
      else if(q < 2.0){
        vec_value = -alpha*((2.0-q)*(2.0-q)/(3.0*q))*dr;
      }
      else{
        vec_value = 0.0;
      }
    }
    return vec_value;
  }
};
