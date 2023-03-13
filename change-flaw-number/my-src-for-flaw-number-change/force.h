#pragma once

#include "prototype.h"

typedef struct{
  PS::F64 p_riemann;
  PS::F64 v_riemann;
  PS::F64 Vij2_hi;
  PS::F64 Vij2_hj;
  PS::F64 ss;
  PS::F64vec dr;
  PS::F64 delta;
  PS::F64vec eij;
  PS::F64vec gradW_hi;
  PS::F64vec gradW_hj;
  PS::matrix Sijast;
} inter_t;

class CalcDensity{
  kernel_t kernel;
 public:
  void operator () (const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens){
    //If we calculate density by summation
    for(PS::S32 i = 0 ; i < Nip ; ++ i){
      dens[i].clear();
      PS::F64 dens_for_smth=0.0;
      for(PS::S32 j = 0 ; j < Njp ; ++ j){
        const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
        dens[i].dens += ep_j[j].mass * kernel.W(dr, ep_i[i].smth);
        dens_for_smth += ep_j[j].mass * kernel.W(dr, PARAM::Csmooth * ep_i[i].smth);
      }
      if(PARAM::IS_VARIABLE_H){
        dens[i].smth = PARAM::SMTH * pow(ep_i[i].mass / dens_for_smth, 1.0/(PS::F64)(PARAM::Dim));
        if(dens[i].smth > PARAM::MAXIMUM_H){
          dens[i].smth = PARAM::MAXIMUM_H;
        }
      }
      else{
        dens[i].smth = ep_i[i].smth;
      }
    }
  }
};

class CalcGradients{
  kernel_t kernel;
 public:
  PS::matrix calc_inverse_tensor_for_3D(const PS::matrix Lab_o){
    const PS::matrix deltaab = PS::matrix(1, 0, 0,
	          0, 1, 0,
	          0, 0, 1 );//Kronecker delta
    PS::matrix Lab;
    PS::F64 det;
    
    det = Lab_o.xx*Lab_o.yy*Lab_o.zz + Lab_o.yx*Lab_o.zy*Lab_o.xz + Lab_o.xy*Lab_o.yz*Lab_o.zx - Lab_o.xz*Lab_o.yy*Lab_o.zx - Lab_o.xy*Lab_o.yx*Lab_o.zz - Lab_o.xx*Lab_o.yz*Lab_o.zy;
    if(fabs(det) < 0.01){
      Lab = deltaab;
    }
    else{
      Lab.xx = (Lab_o.yy*Lab_o.zz-Lab_o.yz*Lab_o.zy)/det;
      Lab.xy = (Lab_o.xz*Lab_o.zy-Lab_o.xy*Lab_o.zz)/det;
      Lab.xz = (Lab_o.xy*Lab_o.yz-Lab_o.xz*Lab_o.yy)/det;
      Lab.yx = (Lab_o.yz*Lab_o.zx-Lab_o.yx*Lab_o.zz)/det;
      Lab.yy = (Lab_o.xx*Lab_o.zz-Lab_o.xz*Lab_o.zx)/det;
      Lab.yz = (Lab_o.xz*Lab_o.yx-Lab_o.xx*Lab_o.yz)/det;
      Lab.zx = (Lab_o.yx*Lab_o.zy-Lab_o.yy*Lab_o.zx)/det;
      Lab.zy = (Lab_o.xy*Lab_o.zx-Lab_o.xx*Lab_o.zy)/det;
      Lab.zz = (Lab_o.xx*Lab_o.yy-Lab_o.xy*Lab_o.yx)/det;
    }

    return(Lab);
  }

  void operator () (const EPI::Gradients* const ep_i, const PS::S32 Nip, const EPJ::Gradients* const ep_j, const PS::S32 Njp, RESULT::Gradients* const gradients){
    const PS::matrix deltaab = PS::matrix(1, 0, 0,
	          0, 1, 0,
	          0, 0, 1 );//Kronecker delta
    for(PS::S32 i = 0 ; i < Nip ; ++ i){
      gradients[i].clear();
      const PS::F64 rhoi_inv2 = 1.0 / ( ep_i[i].dens * ep_i[i].dens );
      for(PS::S32 j = 0 ; j < Njp ; ++ j){
        const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
        const PS::F64vec gradW = kernel.gradW(dr, ep_i[i].smth);
        const PS::F64vec dv = ep_j[j].vel - ep_i[i].vel;
        gradients[i].gradV += -rhoi_inv2 * ep_j[j].mass * gradW;
        gradients[i].graddens += ep_j[j].mass * ( ( ep_j[j].dens - ep_i[i].dens ) / ep_j[j].dens ) * gradW;
        gradients[i].gradpres += ep_j[j].mass * ( ( ep_j[j].pres - ep_i[i].pres ) / ep_j[j].dens ) * gradW;
        gradients[i].gradvel += ( ep_j[j].mass / ep_j[j].dens ) * CalcTensorProduct( dv, gradW );
        if(PARAM::IS_CORRECTED_VELOCITY_GRADIENT==1){
          gradients[i].Lab += ( ep_j[j].mass / ep_j[j].dens ) * CalcTensorProduct( gradW, -dr ); 
        }
      }
      //take inverse of Lab
      if(PARAM::IS_CORRECTED_VELOCITY_GRADIENT==1){
        gradients[i].Lab = calc_inverse_tensor_for_3D(gradients[i].Lab);
      }
      else{
        gradients[i].Lab = deltaab;
      }
    }
  }
};

class CalcElasticForceGodunov{
  const kernel_t kernel;
 public:
  void calc_Vij2_and_ss(const EPI::Elastic ep_i, const EPJ::Elastic ep_j, inter_t* const inter_jP, PS::F64* const ssP, const PS::F64 delta, const PS::F64vec eij){
    const PS::F64 dVi = ep_i.gradV * eij;
    const PS::F64 dVj = ep_j.gradV * eij;
    const PS::F64 Vi = 1.0 / ep_i.dens;
    const PS::F64 Vj = 1.0 / ep_j.dens;
    const PS::F64 hi2 = ep_i.smth * ep_i.smth;
    const PS::F64 hj2 = ep_j.smth * ep_j.smth;
    const PS::F64 rhoij = 0.5 * ( ep_i.dens + ep_j.dens );
    const PS::F64 Vij2crit = 10.0 * ( 1.0 / ( rhoij * rhoij ) );
    const PS::F64 Pij = 0.5 * ( ep_i.pres + ep_j.pres );
    PS::F64 Aij, Bij, Cij, Dij;
    
    //calculate Vij2. According to Sugiura and Inutsuka 2002, in the case of 3 dimension, if Pij is positive we use linear interpolation and if Pij is negative we use cubic spline interpolation to suppress the Tensile Instability.
    if( Pij > 0.0 ){
      //linear interpolation
      Cij = ( Vi - Vj ) / delta;
      Dij = 0.5 * ( Vi + Vj );
      inter_jP->Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
      inter_jP->Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
      *ssP = 0.5 * ( hi2 * Cij * Dij / ( 2 * inter_jP->Vij2_hi ) + hj2 * Cij * Dij / ( 2 * inter_jP->Vij2_hj ) );
      inter_jP->ss = *ssP;
    }
    else{
      //cubic spline interpolation
      Aij = -2 * ( Vi - Vj ) * pow(delta,-3) + ( dVi + dVj ) * pow(delta,-2);
      Bij = 0.5 * ( dVi - dVj ) / delta;
      Cij = 1.5 *( Vi - Vj ) / delta - 0.25 * ( dVi + dVj );
      Dij = 0.5 * ( Vi + Vj ) - 0.125 * ( dVi - dVj ) * delta;
      inter_jP->Vij2_hi = ( ( 0.234375 * hi2 * Aij * Aij + 0.1875 * ( 2 * Aij * Cij + Bij * Bij ) ) * hi2 + 0.25 * ( 2 * Bij * Dij + Cij * Cij ) ) * hi2 + Dij * Dij;
      inter_jP->Vij2_hj = ( ( 0.234375 * hj2 * Aij * Aij + 0.1875 * ( 2 * Aij * Cij + Bij * Bij ) ) * hj2 + 0.25 * ( 2 * Bij * Dij + Cij * Cij ) ) * hj2 + Dij * Dij;
      *ssP = 0.5 * ( ( ( 0.46875 * hi2 * Aij * Bij + 0.375 * ( Aij * Dij + Bij * Cij ) ) * hi2 + 0.5 * Cij * Dij ) * hi2 / inter_jP->Vij2_hi + ( ( 0.46875 * hj2 * Aij * Bij + 0.375 * ( Aij * Dij + Bij * Cij ) ) * hj2 + 0.5 * Cij * Dij ) * hj2 / inter_jP->Vij2_hj );
      inter_jP->ss = *ssP;
    }
    
    //if Vij is too large because of interpolation, we use linear interpolation.
    if( inter_jP->Vij2_hi > Vij2crit || inter_jP->Vij2_hj > Vij2crit ){
      Cij = ( Vi - Vj ) / delta;
      Dij = 0.5 * ( Vi + Vj );
      inter_jP->Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
      inter_jP->Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
      *ssP = 0.5 * ( hi2 * Cij * Dij / ( 2 * inter_jP->Vij2_hi ) + hj2 * Cij * Dij / ( 2 * inter_jP->Vij2_hj ) );
      inter_jP->ss = *ssP;
    }
  }

  //decide riemann solver criterion. If we use riemann solver for elastic equation of state, return 1. If we use that for ideal gas equation of state, return 0.
  int riemann_solver_criterion(const EPI::Elastic ep_i,const EPJ::Elastic ep_j){
    if(PARAM::EQUATION_OF_STATE[ep_i.property_tag]==1||PARAM::EQUATION_OF_STATE[ep_j.property_tag]==1){
      //for elastic equation of state
      return(1);
    }
    else if(PARAM::EQUATION_OF_STATE[ep_i.property_tag]==0||PARAM::EQUATION_OF_STATE[ep_j.property_tag]==0){
      //for ideal gas equation of state
      return(0);
    }
    else{
      //for other equation of state
      if( ep_i.pres<0.0 || ep_j.pres<0.0 || 0.5*(PARAM::U_CRIT_RIEMANN[ep_i.property_tag]+PARAM::U_CRIT_RIEMANN[ep_j.property_tag]) > 0.5*(ep_i.eng+ep_j.eng) ) return(1);
      else return(0);
    }
  }
  
  //solver riemann problem for ideal gas eos
  void solve_riemann_problem_for_ideal_gas_eos(const PS::F64 p1, const PS::F64 rho1, const PS::F64 v1, const PS::F64 gamma1, const PS::F64 p2, const PS::F64 rho2, const PS::F64 v2, const PS::F64 gamma2, inter_t* const inter_jP){
    PS::F64 ppre,p;
    PS::F64 W1,W2;
    const PS::F64 gamma = (gamma1 < gamma2) ? gamma1 : gamma2;//0.5 * ( gamma1 + gamma2 );
    const PS::F64 alpha = ( 2.0 * gamma ) / ( gamma - 1.0 );
    
    p = 0.5 * ( p1 + p2 );
    for( PS::U32 loop=0 ; loop<5 ; loop++){
      ppre = p;
      if( p < p1 - PARAM::ACC ) W1 = sqrt( p1 * rho1 ) * ( ( gamma - 1.0 ) / ( 2.0 * sqrt(gamma) ) ) * ( 1.0 - p / p1 ) / ( 1.0 - pow(p/p1 , 1.0/alpha) );
      else W1 = sqrt(p1 * rho1) * sqrt( 0.5 * ( gamma + 1.0 ) * p / p1 + 0.5 * ( gamma - 1.0 ) );
      if( p < p2 - PARAM::ACC ) W2 = sqrt( p2 * rho2 ) * ( ( gamma - 1.0 ) / ( 2.0 * sqrt(gamma) ) ) * ( 1.0 - p / p2 ) / ( 1.0 - pow(p/p2 , 1.0/alpha) );
      else W2 = sqrt(p2 * rho2) * sqrt( 0.5 * ( gamma + 1.0 ) * p / p2 + 0.5 * ( gamma - 1.0 ) );
      p = ( ( p2 / W2 + p1 / W1 ) + v2 - v1 ) / ( 1.0 / W2 + 1.0 / W1 );
      if( fabs( p - ppre ) < PARAM::ACC ) break;
    }
    if( std::isnan(p) || std::isinf(p) ){
      inter_jP->p_riemann = 0.5 * ( p1 + p2 );
      inter_jP->v_riemann = 0.5 * ( v1 + v2 );
    }
    else{
      inter_jP->p_riemann = p;
      inter_jP->v_riemann = ( ( W1 * v1 + W2 * v2 ) + p2 - p1 ) / ( W1 + W2 );
    }
  }

  //solve riemann problem for elastic EoS of P=Cs^2(rho-rho0eos)
  void solve_riemann_problem_for_elastic_eos(const PS::F64 p1,const PS::F64 rho1,const PS::F64 v1,const PS::F64 Csi,const PS::F64 p2,const PS::F64 rho2,const PS::F64 v2,const PS::F64 Csj,inter_t* const inter_jP)
  {
    PS::F64 ppre,p;
    PS::F64 W1,W2;	
    PS::F64 Csij = 0.5*(Csi+Csj);
    PS::F64 rho0eosij = 0.5*( ( rho1 - p1/(Csij*Csij) ) + ( rho2 - p2/(Csij*Csij) ) );
    
    p = 0.5 * ( p1 + p2 );
    for(PS::S32 loop=0 ; loop<5 ; loop++){
      ppre = p;
      if( p < p1 - PARAM::ACC ) W1 = fabs( ( (p-p1) / Csij ) * ( 1.0 / log( Csij*Csij*rho1 / (p+Csij*Csij*rho0eosij) ) ) );
      else W1 = Csij * sqrt( rho1 * ( rho1 + (p-p1)/(Csij*Csij) ) );
      if( p < p2 - PARAM::ACC ) W2 = fabs( ( (p-p2) / Csij ) * ( 1.0 / log( Csij*Csij*rho2 / (p+Csij*Csij*rho0eosij) ) ) );
      else W2 = Csij * sqrt( rho2 * ( rho2 + (p-p2)/(Csij*Csij) ) );
      p = ( ( p2 / W2 + p1 / W1 ) + v2 - v1 ) / ( 1.0 / W2 + 1.0 / W1 );
      if( fabs( p - ppre ) < PARAM::ACC ) break;
    }
    if( std::isnan(p) || std::isinf(p) ){
      inter_jP->p_riemann = 0.5 * ( p1 + p2 );
      inter_jP->v_riemann = 0.5 * ( v1 + v2 );
    }
    else{
      inter_jP->p_riemann = p;
      inter_jP->v_riemann = 0.5 * ( v1 + v2 );
    }
  }
  
  void calc_riemann_solver(const EPI::Elastic ep_i, const EPJ::Elastic ep_j, inter_t* const inter_jP, const PS::F64 ss, const PS::F64 delta, const PS::F64vec eij, const PS::F64 dt){
    PS::F64 drdsi=0.0,drdsj=0.0,dpdsi=0.0,dpdsj=0.0,dvdsi=0.0,dvdsj=0.0,vi,vj;
    PS::F64 rhoR,rhoL,pR,pL,vR,vL;
    
    vi = ep_i.vel * eij;
    vj = ep_j.vel * eij;
    
    //if space order is second or particles are not in shock front, calculate s direction's gradients 
    if( ( PARAM::SPACE_ORDER==2 ) && ( 3*(vj-vi) < ( ( ep_i.snds < ep_j.snds ) ? ep_i.snds : ep_j.snds ) ) ){
      drdsi = ep_i.graddens * eij;
      drdsj = ep_j.graddens * eij;
      dpdsi = ep_i.gradpres * eij;
      dpdsj = ep_j.gradpres * eij;
      dvdsi = CalcTensorToScalar( ep_i.gradvel, CalcTensorProduct( eij, eij ) );
      dvdsj = CalcTensorToScalar( ep_j.gradvel, CalcTensorProduct( eij, eij ) );
      //monotonisity constraint
      if( dvdsi * dvdsj < 0.0 ) dvdsi = dvdsj = 0.0;
      if( drdsi * drdsj < 0.0 ) drdsi = drdsj = 0.0;
      if( dpdsi * dpdsj < 0.0 ) dpdsi = dpdsj = 0.0;
    }
    
    //set initial condision of riemann problem
    rhoR = ep_i.dens + drdsi * ( ss + ep_i.snds * 0.5 * dt - 0.5 * delta );
    pR = ep_i.pres + dpdsi * ( ss + ep_i.snds * 0.5 * dt - 0.5 * delta );
    vR = vi + dvdsi * ( ss + ep_i.snds * 0.5 * dt - 0.5 * delta );
    rhoL = ep_j.dens + drdsj * ( ss - ep_j.snds * 0.5 * dt + 0.5 * delta );
    pL = ep_j.pres + dpdsj * ( ss - ep_j.snds * 0.5 * dt + 0.5 * delta );
    vL = vj + dvdsj * ( ss - ep_j.snds * 0.5 * dt + 0.5 * delta );
    if( rhoR < 0.0 ) rhoR = ep_i.dens;
    if( rhoL < 0.0 ) rhoL = ep_j.dens;
    if( pR < 0.0 && PARAM::EQUATION_OF_STATE[ep_i.property_tag]==0 ) pR = ep_i.pres;
    if( pL < 0.0 && PARAM::EQUATION_OF_STATE[ep_j.property_tag]==0 ) pL = ep_j.pres;
    
    //solve riemann solver
    if( riemann_solver_criterion(ep_i, ep_j) ){
      solve_riemann_problem_for_elastic_eos(pR,rhoR,vR,ep_i.snds,pL,rhoL,vL,ep_j.snds,inter_jP);
    }
    else{
      solve_riemann_problem_for_ideal_gas_eos(pR,rhoR,vR,ep_i.gamma,pL,rhoL,vL,ep_j.gamma,inter_jP);
    }
  }

  PS::matrix calc_epsilon_dot_rho_godunov(const PS::F64 massj, const inter_t inter_j, const PS::F64vec v_ast, const PS::F64vec tcenter_v){
    PS::matrix epsilon_dot_rho;
    epsilon_dot_rho.xx = massj * (v_ast.x-tcenter_v.x) * ( inter_j.Vij2_hi * inter_j.gradW_hi.x + inter_j.Vij2_hj * inter_j.gradW_hj.x );
    epsilon_dot_rho.yy = massj * (v_ast.y-tcenter_v.y) * ( inter_j.Vij2_hi * inter_j.gradW_hi.y + inter_j.Vij2_hj * inter_j.gradW_hj.y );
    epsilon_dot_rho.zz = massj * (v_ast.z-tcenter_v.z) * ( inter_j.Vij2_hi * inter_j.gradW_hi.z + inter_j.Vij2_hj * inter_j.gradW_hj.z );
    epsilon_dot_rho.xy = epsilon_dot_rho.yx = 0.5 * massj * ( (v_ast.x-tcenter_v.x) * ( inter_j.Vij2_hi * inter_j.gradW_hi.y + inter_j.Vij2_hj * inter_j.gradW_hj.y ) 
	                              + (v_ast.y-tcenter_v.y) * ( inter_j.Vij2_hi * inter_j.gradW_hi.x + inter_j.Vij2_hj * inter_j.gradW_hj.x ) );
    epsilon_dot_rho.xz = epsilon_dot_rho.zx = 0.5 * massj * ( (v_ast.x-tcenter_v.x) * ( inter_j.Vij2_hi * inter_j.gradW_hi.z + inter_j.Vij2_hj * inter_j.gradW_hj.z ) 
	                              + (v_ast.z-tcenter_v.z) * ( inter_j.Vij2_hi * inter_j.gradW_hi.x + inter_j.Vij2_hj * inter_j.gradW_hj.x ) );
    epsilon_dot_rho.yz = epsilon_dot_rho.zy = 0.5 * massj * ( (v_ast.y-tcenter_v.y) * ( inter_j.Vij2_hi * inter_j.gradW_hi.z + inter_j.Vij2_hj * inter_j.gradW_hj.z ) 
	                              + (v_ast.z-tcenter_v.z) * ( inter_j.Vij2_hi * inter_j.gradW_hi.y + inter_j.Vij2_hj * inter_j.gradW_hj.y ) );
    
    return(epsilon_dot_rho);
  }
  
  PS::matrix calc_R_rho_godunov(const PS::F64 massj, const inter_t inter_j, const PS::F64vec v_ast, const PS::F64vec tcenter_v){
    PS::matrix R_rho;
    R_rho.xx = R_rho.yy = R_rho.zz = 0.0;
    R_rho.xy = 0.5 * massj * ( (v_ast.x-tcenter_v.x) * ( inter_j.Vij2_hi * inter_j.gradW_hi.y + inter_j.Vij2_hj * inter_j.gradW_hj.y ) 
                               - (v_ast.y-tcenter_v.y) * ( inter_j.Vij2_hi * inter_j.gradW_hi.x + inter_j.Vij2_hj * inter_j.gradW_hj.x ) );
    R_rho.yx = -R_rho.xy;
    R_rho.xz = 0.5 * massj * ( (v_ast.x-tcenter_v.x) * ( inter_j.Vij2_hi * inter_j.gradW_hi.z + inter_j.Vij2_hj * inter_j.gradW_hj.z ) 
                               - (v_ast.z-tcenter_v.z) * ( inter_j.Vij2_hi * inter_j.gradW_hi.x + inter_j.Vij2_hj * inter_j.gradW_hj.x ) );
    R_rho.zx = -R_rho.xz;
    R_rho.yz = 0.5 * massj * ( (v_ast.y-tcenter_v.y) * ( inter_j.Vij2_hi * inter_j.gradW_hi.z + inter_j.Vij2_hj * inter_j.gradW_hj.z ) 
                               - (v_ast.z-tcenter_v.z) * ( inter_j.Vij2_hi * inter_j.gradW_hi.y + inter_j.Vij2_hj * inter_j.gradW_hj.y ) );
    R_rho.zy = -R_rho.yz;
    
    return(R_rho);
  }
  
  void operator () (const EPI::Elastic* const ep_i, const PS::S32 Nip, const EPJ::Elastic* const ep_j, const PS::S32 Njp, RESULT::Elastic* const elast){
    const PS::F64 dt = getTimeStep();
    const PS::matrix deltaab = PS::matrix(1, 0, 0,
	          0, 1, 0,
	          0, 0, 1 );//Kronecker delta
    inter_t* inter;
    inter = (inter_t*)malloc( Njp * sizeof(inter_t) );
    PS::F64vec tcenter_v;
    for(PS::S32 i = 0; i < Nip ; ++ i){
      elast[i].clear();
      for(PS::S32 j = 0 ; j < Njp ; ++ j){
        if( ep_i[i].id != ep_j[j].id ){
          inter[j].dr = ep_i[i].pos - ep_j[j].pos;
          inter[j].delta = sqrt( inter[j].dr * inter[j].dr );
          inter[j].eij = inter[j].dr / inter[j].delta;
          PS::F64 ss;
          
          //calculate Vij2
          calc_Vij2_and_ss(ep_i[i],ep_j[j],&(inter[j]),&ss,inter[j].delta,inter[j].eij);
          inter[j].Sijast = 0.5 * (ep_i[i].Sab + ep_j[j].Sab) + (ep_i[i].Sab - ep_j[j].Sab) * (ss/inter[j].delta);
          
          //calculate riemann solver
          calc_riemann_solver(ep_i[i],ep_j[j],&(inter[j]),ss,inter[j].delta,inter[j].eij,dt);
          
          inter[j].gradW_hi = kernel.gradW(inter[j].dr, sqrt(2.0) * ep_i[i].smth);
          inter[j].gradW_hj = kernel.gradW(inter[j].dr, sqrt(2.0) * ep_j[j].smth);
          
          elast[i].acc += ( - ep_j[j].mass * inter[j].p_riemann * ( inter[j].Vij2_hi * inter[j].gradW_hi + inter[j].Vij2_hj * inter[j].gradW_hj ) )
            + ( ep_j[j].mass * inter[j].Sijast * ( inter[j].Vij2_hi * inter[j].gradW_hi + inter[j].Vij2_hj * inter[j].gradW_hj ) );
        }
      }
      tcenter_v = ep_i[i].vel + 0.5 * dt * elast[i].acc;
      for(PS::S32 j = 0 ; j < Njp ; ++ j){
        if( ep_i[i].id != ep_j[j].id ){
          PS::F64vec v_ast;
          if( PARAM::EQUATION_OF_STATE[ep_i[i].property_tag]==0 && PARAM::EQUATION_OF_STATE[ep_j[j].property_tag]==0 ){
            const PS::F64vec vij_vector = 0.5 * ( ep_i[i].vel + ep_j[j].vel ) + ( ep_i[i].vel - ep_j[j].vel ) * ( inter[j].ss / inter[j].delta);
            const PS::F64 vij = vij_vector * inter[j].eij;
            v_ast = inter[j].v_riemann * inter[j].eij + ( vij_vector - vij * inter[j].eij );
          }
          else{
            v_ast = 0.5 * ( ep_i[i].vel + ep_j[j].vel ) + ( ep_i[i].vel - ep_j[j].vel ) * ( inter[j].ss / inter[j].delta);
          }
          const PS::matrix epsilon_dot_rho = calc_epsilon_dot_rho_godunov(ep_j[j].mass, inter[j], v_ast, tcenter_v);
          const PS::matrix R_rho = calc_R_rho_godunov(ep_j[j].mass, inter[j], v_ast, tcenter_v);
          
          elast[i].eng_dot += CalcTensorToScalar( -inter[j].p_riemann * deltaab + inter[j].Sijast, epsilon_dot_rho );
          elast[i].dSab_rho_dt += 2.0*PARAM::MU_SHEAR[ep_i[i].property_tag] * ( epsilon_dot_rho - 0.333333333333 * epsilon_dot_rho.getTrace() * deltaab ) + inter[j].Sijast * R_rho + R_rho * inter[j].Sijast + inter[j].Sijast * epsilon_dot_rho.getTrace();
          if(PARAM::DENSITY_DEVELOPMENT_METHOD==1){
            elast[i].drho_dt += ep_j[j].mass * ( ep_i[i].dens / ep_j[j].dens ) * ( ep_i[i].vel - ep_j[j].vel ) * kernel.gradW(inter[j].dr, ep_i[i].smth);
          }
        }
      }
    }

    free(inter);
  }
};

class CalcElasticForceStand{
  const kernel_t kernel;
 public:
  PS::matrix calc_epsilon_dot_rho_stand(const EPI::Elastic ep_i, const EPJ::Elastic ep_j, const PS::F64vec gradW){
    PS::matrix epsilon_dot_rho;
    epsilon_dot_rho.xx = (ep_j.mass/(ep_i.dens*ep_j.dens)) * (ep_j.vel.x-ep_i.vel.x) * gradW.x;
    epsilon_dot_rho.yy = (ep_j.mass/(ep_i.dens*ep_j.dens)) * (ep_j.vel.y-ep_i.vel.y) * gradW.y;
    epsilon_dot_rho.zz = (ep_j.mass/(ep_i.dens*ep_j.dens)) * (ep_j.vel.z-ep_i.vel.z) * gradW.z;
    epsilon_dot_rho.xy = epsilon_dot_rho.yx = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.x-ep_i.vel.x) * gradW.y + (ep_j.vel.y-ep_i.vel.y) * gradW.x );
    epsilon_dot_rho.xz = epsilon_dot_rho.zx = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.x-ep_i.vel.x) * gradW.z + (ep_j.vel.z-ep_i.vel.z) * gradW.x );
    epsilon_dot_rho.yz = epsilon_dot_rho.zy = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.y-ep_i.vel.y) * gradW.z + (ep_j.vel.z-ep_i.vel.z) * gradW.y );

    return(epsilon_dot_rho);
  }

  PS::matrix calc_R_rho_stand(const EPI::Elastic ep_i, const EPJ::Elastic ep_j, const PS::F64vec gradW){
    PS::matrix R_rho;
    R_rho.xx = R_rho.yy = R_rho.zz = 0.0;
    R_rho.xy = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.x-ep_i.vel.x) * gradW.y - (ep_j.vel.y-ep_i.vel.y) * gradW.x );
    R_rho.yx = -R_rho.xy;
    R_rho.xz = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.x-ep_i.vel.x) * gradW.z - (ep_j.vel.z-ep_i.vel.z) * gradW.x );
    R_rho.zx = -R_rho.xz;
    R_rho.yz = (ep_j.mass/(2.0*ep_i.dens*ep_j.dens)) * ( (ep_j.vel.y-ep_i.vel.y) * gradW.z - (ep_j.vel.z-ep_i.vel.z) * gradW.y );
    R_rho.zy = -R_rho.yz;

    return(R_rho);
  }

  void operator () (const EPI::Elastic* const ep_i, const PS::S32 Nip, const EPJ::Elastic* const ep_j, const PS::S32 Njp, RESULT::Elastic* const elast){
    for(PS::S32 i = 0; i < Nip ; ++ i){
      elast[i].clear();
      const PS::matrix deltaab = PS::matrix(1, 0, 0,
	            0, 1, 0,
	            0, 0, 1 );//Kronecker delta
      for(PS::S32 j = 0; j < Njp ; ++ j){
        const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
        const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
        const PS::F64 rho_ij = 0.5 * ( ep_i[i].dens + ep_j[j].dens );
        const PS::F64 Cs_ij = 0.5 * ( ep_i[i].snds + ep_j[j].snds );
        const PS::F64 h_ij = 0.5 * ( ep_i[i].smth + ep_j[j].smth );
        PS::F64 mu_ij = ( h_ij * dv * dr ) / ( dr * dr + 0.01 * h_ij * h_ij );
        if(mu_ij>0.0) mu_ij=0.0;
        const PS::F64 AV = - PARAM::ALPHA_VIS * ( Cs_ij * mu_ij / rho_ij ) + PARAM::BETA_VIS * ( mu_ij * mu_ij / rho_ij );
        const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ep_i[i].smth) + kernel.gradW(dr, ep_j[j].smth));
        PS::F64vec gradW_hi;
        if(PARAM::IS_CORRECTED_VELOCITY_GRADIENT==1){
          gradW_hi = ep_i[i].Lab * kernel.gradW(dr, ep_i[i].smth);
        }
        else{
          gradW_hi = kernel.gradW(dr, ep_i[i].smth);
        }
        PS::F64vec gradW_hi_rho = kernel.gradW(dr, ep_i[i].smth);
        const PS::matrix epsilon_dot_rho = calc_epsilon_dot_rho_stand(ep_i[i], ep_j[j], gradW_hi);
        const PS::matrix R_rho = calc_R_rho_stand(ep_i[i], ep_j[j], gradW_hi);
        const PS::matrix epsilon_dot_rho_eng = calc_epsilon_dot_rho_stand(ep_i[i], ep_j[j], gradW);
        elast[i].acc     += ( -ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * gradW)
          + ( ep_j[j].mass * (ep_i[i].Sab / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].Sab / (ep_j[j].dens * ep_j[j].dens) ) * gradW);
        elast[i].eng_dot += (0.5 * ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * dv * gradW)
          + CalcTensorToScalar(ep_i[i].Sab, epsilon_dot_rho_eng);
        elast[i].dSab_rho_dt += 2.0*PARAM::MU_SHEAR[ep_i[i].property_tag] * ( epsilon_dot_rho - 0.333333333333 * epsilon_dot_rho.getTrace() * deltaab ) + ep_i[i].Sab * R_rho + R_rho * ep_i[i].Sab + ep_i[i].Sab * epsilon_dot_rho.getTrace();
        if(PARAM::DENSITY_DEVELOPMENT_METHOD==1){
          elast[i].drho_dt += ep_j[j].mass * ( ep_i[i].dens / ep_j[j].dens ) * dv * gradW_hi_rho;
        }
      }
    }
  }
};

class CalcGravityForceCloseParticle{
 public:
  void operator () (const EPI::Grav* const ep_i, const PS::S32 Nip, const EPJ::Grav* const ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
    for(PS::S32 i = 0; i < Nip ; ++ i){
      grav[i].clear();
      for(PS::S32 j = 0; j < Njp ; ++ j){
        if(ep_i[i].id != ep_j[j].id){
          const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
          const PS::F64 dr3_inv = pow(dr*dr,-1.5);
          const PS::F64 delta = sqrt(dr*dr);
          const PS::F64 hmax = ( ( ep_i[i].smth > ep_j[j].smth ) ? ep_i[i].smth : ep_j[j].smth );
          //if two particles are far away
          if( delta > 2.0*hmax ){ 
            grav[i].grav += - PARAM::G_CONST * ep_j[j].mass * dr * dr3_inv;
          }   
          else{
            const PS::F64 hratio = ( ( ep_i[i].smth < ep_j[j].smth ) ? ep_i[i].smth : ep_j[j].smth ) / hmax;
            const PS::F64 Xtilda = ( delta / hmax ) * (1.0 / (1.0 + hratio * (0.007697 + hratio * (0.53799 - hratio * 0.125256))));
            grav[i].grav += - PARAM::G_CONST * ep_j[j].mass * dr * dr3_inv * ( erf(Xtilda) - (2.0 / sqrt(M_PI)) * Xtilda * exp(- Xtilda * Xtilda));
          }
        }
      }
    }
  }
};

class CalcGravityForceSuperParticle{
 public:
  void operator () (const EPI::Grav* const ep_i, const PS::S32 Nip, const PS::SPJMonopole* const sp_j, const PS::S32 Njp, RESULT::Grav* const grav){
    for(PS::S32 i = 0; i < Nip ; ++ i){
      for(PS::S32 j = 0; j < Njp ; ++ j){
        const PS::F64vec dr = ep_i[i].pos - sp_j[j].pos;
        const PS::F64 dr3_inv = pow(dr*dr,-1.5);
        grav[i].grav += - PARAM::G_CONST * sp_j[j].mass * dr * dr3_inv;
      }
    }
  }
};
