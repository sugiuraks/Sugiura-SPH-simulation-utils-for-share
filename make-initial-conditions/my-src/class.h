#pragma once

class FileHeader{
 public:
  int Nbody;
  double time;
  int max_flaw_number;
  int readAscii(FILE* fp){
    fscanf(fp, "%le\n", &time);
    fscanf(fp, "%d\n", &Nbody);
    return Nbody;
  }
  void writeAscii(FILE* fp) const{
    fprintf(fp, "%e\n", time);
    fprintf(fp, "%d\n", Nbody);
  }
  int readBinary(FILE* fp){
    fread(&time, sizeof(double), 1, fp);
    fread(&Nbody, sizeof(int), 1, fp);
    fread(&max_flaw_number, sizeof(int), 1, fp);
    if(max_flaw_number!=MAX_FLAW_NUMBER){
      printf("max_flaw_number of this input file is different from defined MAX_FLAW_NUMBER! \n");
      printf("max_flaw_number of this input file: %d defined MAX_FLAW_NUMBER: %d \n",max_flaw_number,MAX_FLAW_NUMBER);
    }
    return Nbody;
  }
  void writeBinary(FILE* fp) const{
    fwrite(&time, sizeof(double), 1, fp);
    fwrite(&Nbody, sizeof(int), 1, fp);
    int temp=MAX_FLAW_NUMBER;
    fwrite(&temp, sizeof(int), 1, fp);
  }
};

namespace RESULT{
  //Density summation
  class Dens{
  public:
    PS::F64 dens;
    PS::F64 smth;
    void clear(){
      dens = 0;
    }
  };
  class Gradients{
  public:
    PS::F64vec gradV;
    PS::F64vec graddens;
    PS::F64vec gradpres;
    PS::matrix gradvel;
    PS::matrix Lab;
    void clear(){
      gradV = 0.0;
      graddens = 0.0;
      gradpres = 0.0;
      gradvel = 0.0;
      Lab = 0.0;
    }
  };  
  class Elastic{
  public:
    PS::F64vec acc;
    PS::F64 eng_dot;
    PS::matrix dSab_rho_dt;
    PS::F64 drho_dt;
    void clear(){
      acc = 0;
      eng_dot = 0;
      dSab_rho_dt = 0;
      drho_dt = 0;
    }
  };
  //Self gravity
  class Grav{
  public:
    PS::F64vec grav;
    void clear(){
      grav = 0.0;
    }
  };
}

class RealPtcl{
 public:
  PS::F64 mass; //mass of particle 
  PS::F64vec pos, vel, acc; //VELocity, POSition, ACCeleration 
  PS::F64vec grav; //self GRAVity
  PS::F64vec acc_old; //acceleration of previous time step
  PS::F64vec grav_old; //self gravity of previous time step
  PS::F64 dens;//DENSity
  PS::F64 eng; //ENerGy
  PS::F64 pres;//PRESsure
  PS::F64 smth;//SMooTHing length
  PS::F64 snds; //SouND Speed
  PS::F64 temp; //TEMPerature

  PS::F64 eng_max; //max internal energy that this particle experienced
  
  PS::F64 eng_dot; //time change rate of internal energy
  PS::F64 eng_dot_old; //time change rate of internal energy of previous time step
  PS::F64vec vel_old; //velocity of previous time step
  PS::F64vec pos_old; //position of previous time step
  PS::F64 eng_old; //internal energy of previous time step   
  PS::F64 dt; //time step for this particle
  PS::S64 id; //ID
  PS::S32 property_tag; //tag representing the type of material for this particle
  PS::F64 gamma; //effective gamma for riemann solver of ideal gas EoS
  
  PS::F64vec gradV; //gradient of specific volume
  PS::F64vec graddens; //gradient of density
  PS::F64vec gradpres; //gradient of pressure
  PS::matrix gradvel; //gradient matrix of velocity

  PS::matrix Sab; //deviatoric stress tensor
  PS::matrix Sab_rho; //deviatoric stress tensor / density
  PS::matrix Sab_rho_old; //deviatoric stress tensor  / density of previous time step
  PS::matrix dSab_rho_dt; //time derivative of Sab/rho
  PS::matrix dSab_rho_dt_old; //time derivative of Sab/rho of previous time step

  PS::F64 drho_dt; //time derivative of density
  PS::F64 drho_dt_old; //time derivative of density of previous time step
  PS::F64 dens_old; //density of previous time step

  PS::F64 damage; //damage parameter D
  PS::F64 D_cbrt; //cubic root of damage parameter D
  PS::F64 dD_cbrt_dt; //time change rate of cubic root of D
  PS::F64 eps_ij_act_flaw[MAX_FLAW_NUMBER]; //array of flaw activation threshold strain
  PS::S64 ni_tot_flaw; //the number of flaws that are assigned to this particle
  PS::F64 eps_i_eff; //effective strain 

  PS::F64 alpha_por; //distension parameter
  PS::F64 dalpha_dt; //time change rate of distension parameter
  PS::F64 pres_old; //pressure of previous time step for porosity model

  PS::F64 f_test;
  PS::matrix Lab; //correction term for velocity gradient
  //Copy functions
  void copyFromForce(const RESULT::Dens& dens){
    this->dens = dens.dens;
    this->smth = dens.smth;
  }
  void copyFromForce(const RESULT::Gradients& gradients){
    this->gradV = gradients.gradV;
    this->graddens = gradients.graddens;
    this->gradpres = gradients.gradpres;
    this->gradvel = gradients.gradvel;
    this->Lab = gradients.Lab;
  }
  void copyFromForce(const RESULT::Elastic& force){
    this->acc     = force.acc;
    this->eng_dot = force.eng_dot;
    this->dSab_rho_dt = force.dSab_rho_dt;
    this->drho_dt = force.drho_dt;
  }
  void copyFromForce(const RESULT::Grav& force){
    this->grav = force.grav;
  }
  //Give necessary values to FDPS
  PS::F64 getCharge() const{
    return this->mass;
  }
  PS::F64vec getPos() const{
    return this->pos;
  }
  PS::F64 getRSearch() const{
    PS::F64 value;
    if(PARAM::KERNEL_FUNCTION==0){
      //gaussian kernel
      value = PARAM::SEARCH_RADIUS_H;
    }
    if(PARAM::KERNEL_FUNCTION==1){
      //cubic spline kernel
      value = 2.0;
    }
    return value * this->smth;
  }
  void setPos(const PS::F64vec& pos){
    this->pos = pos;
  }
  void writeAscii(FILE* fp) const{
    double km=1.0e5;
    //fprintf(fp, "%f %f %f %f %f %f %f \n", this->pos.x, this->pos.y, this->pos.z, this->dens, this->pres, this->Sab.xy, this->Sab.yz);
    //fprintf(fp,"%d %e %e %e %e %e %e %e \n",this->id,this->pos.x,this->vel.x,this->eng,this->Sab_rho.xy,this->Sab.xy,this->dens,this->snds);
    fprintf(fp,"%.2f %.2f %.2f %.2f %lld %d \n",this->pos.x/km,this->pos.y/km,this->pos.z/km,this->smth/km,this->id,this->property_tag);
  }
  void readAscii(FILE* fp){
    fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->dens, &this->eng, &this->pres);
  }
  void writeBinary(FILE* fp) const{
    fwrite(&this->pos.x        , sizeof(PS::F64), 1, fp);
    fwrite(&this->pos.y        , sizeof(PS::F64), 1, fp);
    fwrite(&this->pos.z        , sizeof(PS::F64), 1, fp);
    fwrite(&this->vel.x        , sizeof(PS::F64), 1, fp);
    fwrite(&this->vel.y        , sizeof(PS::F64), 1, fp);
    fwrite(&this->vel.z        , sizeof(PS::F64), 1, fp);
    fwrite(&this->mass         , sizeof(PS::F64), 1, fp);
    fwrite(&this->eng          , sizeof(PS::F64), 1, fp);
    fwrite(&this->smth         , sizeof(PS::F64), 1, fp);
    fwrite(&this->dens         , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.xx   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.xy   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.xz   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.yx   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.yy   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.yz   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.zx   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.zy   , sizeof(PS::F64), 1, fp);
    fwrite(&this->Sab_rho.zz   , sizeof(PS::F64), 1, fp);
    fwrite(&this->property_tag , sizeof(PS::S32), 1, fp);
    fwrite(&this->id           , sizeof(PS::S32), 1, fp);
    fwrite(&this->D_cbrt       , sizeof(PS::F64), 1, fp);
    fwrite(&this->alpha_por    , sizeof(PS::F64), 1, fp);
    fwrite(&this->ni_tot_flaw  , sizeof(PS::S64), 1, fp);
    fwrite(this->eps_ij_act_flaw,sizeof(PS::F64),MAX_FLAW_NUMBER,fp);
	fwrite(&this->eng_max      , sizeof(PS::F64), 1, fp);
  }
  void readBinary(FILE* fp){
    fread(&this->pos.x        , sizeof(PS::F64), 1, fp);
    fread(&this->pos.y        , sizeof(PS::F64), 1, fp);
    fread(&this->pos.z        , sizeof(PS::F64), 1, fp);
    fread(&this->vel.x        , sizeof(PS::F64), 1, fp);
    fread(&this->vel.y        , sizeof(PS::F64), 1, fp);
    fread(&this->vel.z        , sizeof(PS::F64), 1, fp);
    fread(&this->mass         , sizeof(PS::F64), 1, fp);
    fread(&this->eng          , sizeof(PS::F64), 1, fp);
    fread(&this->smth         , sizeof(PS::F64), 1, fp);
    fread(&this->dens         , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.xx   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.xy   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.xz   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.yx   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.yy   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.yz   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.zx   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.zy   , sizeof(PS::F64), 1, fp);
    fread(&this->Sab_rho.zz   , sizeof(PS::F64), 1, fp);
    fread(&this->property_tag , sizeof(PS::S32), 1, fp);
    fread(&this->id           , sizeof(PS::S32), 1, fp);
    fread(&this->D_cbrt       , sizeof(PS::F64), 1, fp);
    this->damage = pow(this->D_cbrt, 3.0);
    fread(&this->alpha_por    , sizeof(PS::F64), 1, fp);
    fread(&this->ni_tot_flaw  , sizeof(PS::S64), 1, fp);
    fread(this->eps_ij_act_flaw,sizeof(PS::F64),MAX_FLAW_NUMBER,fp);
	fread(&this->eng_max      , sizeof(PS::F64), 1, fp);
  }
};

namespace EPI{
  class Dens{
  public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64    smth;
    void copyFromFP(const RealPtcl& rp){
      this->pos  = rp.pos;
      this->mass = rp.mass;
      this->smth = rp.smth;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
    void setPos(const PS::F64vec& pos){
      this->pos = pos;
    }
  };
  class Gradients{
  public:
    PS::F64vec pos;
    PS::F64    dens;
    PS::F64    smth;
    PS::F64    pres;
    PS::F64vec vel;
    void copyFromFP(const RealPtcl& rp){
      this->pos  = rp.pos;
      this->dens = rp.dens;
      this->smth = rp.smth;
      this->pres = rp.pres;
      this->vel  = rp.vel;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
    void setPos(const PS::F64vec& pos){
      this->pos = pos;
    }
  };
  class Elastic{
  public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec vel_old;
    PS::F64    smth;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    snds;
    PS::F64    eng;
    PS::S64    id;///DEBUG
    PS::F64vec gradV;
    PS::F64vec graddens;
    PS::F64vec gradpres;
    PS::matrix gradvel;
    PS::S32 property_tag; 
    PS::F64 gamma; 
    PS::matrix Sab;
    PS::matrix Lab;
    void copyFromFP(const RealPtcl& rp){
      this->pos   = rp.pos;
      this->vel   = rp.vel;
      this->vel_old = rp.vel_old;
      this->smth  = rp.smth;
      this->dens  = rp.dens;
      this->pres  = rp.pres;
      this->snds = rp.snds;
      this->eng = rp.eng;
      this->id   = rp.id;///DEBUG
      this->gradV = rp.gradV;
      this->graddens = rp.graddens;
      this->gradpres = rp.gradpres;
      this->gradvel = rp.gradvel;
      this->property_tag = rp.property_tag;
      this->gamma = rp.gamma;
      this->Sab = rp.Sab;
      this->Lab = rp.Lab;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::IS_GODUNOV_SPH){
        //Godunov SPH with gaussian kernel
        return sqrt( 2.0 ) * PARAM::SEARCH_RADIUS_H * this->smth;
      }
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
  };
  class Grav{
  public:
    PS::F64vec pos;
    PS::F64    smth;
    PS::S64    id;
    PS::F64vec getPos() const{
      return this->pos;
    }
    void copyFromFP(const RealPtcl& rp){
      this->pos = rp.pos;
      this->smth = rp.smth;
      this->id  = rp.id;
    }
  };
}

namespace EPJ{
  class Dens{
  public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64    smth;
    void copyFromFP(const RealPtcl& rp){
      this->mass = rp.mass;
      this->pos  = rp.pos;
      this->smth = rp.smth;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
    void setPos(const PS::F64vec& pos){
      this->pos = pos;
    }
  };
  class Gradients{
  public:
    PS::F64vec pos;
    PS::F64    dens;
    PS::F64vec vel;
    PS::F64    pres;
    PS::F64    mass;
    PS::F64    smth;
    void copyFromFP(const RealPtcl& rp){
      this->pos  = rp.pos;
      this->dens = rp.dens;
      this->vel  = rp.vel;
      this->pres = rp.pres;
      this->mass = rp.mass;
      this->smth = rp.smth;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
    void setPos(const PS::F64vec& pos){
      this->pos = pos;
    }
  };
  class Elastic{
  public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    smth;
    PS::F64    mass;
    PS::F64    snds;
    PS::F64    eng;
    PS::S64    id;///DEBUG
    PS::F64vec gradV;
    PS::F64vec graddens;
    PS::F64vec gradpres;
    PS::matrix gradvel;
    PS::S32 property_tag; 
    PS::F64 gamma;
    PS::matrix Sab; 
    void copyFromFP(const RealPtcl& rp){
      this->pos   = rp.pos;
      this->vel   = rp.vel;
      this->dens  = rp.dens;
      this->pres  = rp.pres;
      this->smth  = rp.smth;
      this->mass  = rp.mass;
      this->snds  = rp.snds;
      this->eng = rp.eng;
      this->id    = rp.id;///DEBUG
      this->gradV = rp.gradV;
      this->graddens = rp.graddens;
      this->gradpres = rp.gradpres;
      this->gradvel = rp.gradvel;
      this->property_tag = rp.property_tag;
      this->gamma = rp.gamma;
      this->Sab = rp.Sab;
    }
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getRSearch() const{
      PS::F64 value;
      if(PARAM::IS_GODUNOV_SPH){
        //Godunov SPH with gaussian kernel
        return sqrt( 2.0 ) * PARAM::SEARCH_RADIUS_H * this->smth;
      }
      if(PARAM::KERNEL_FUNCTION==0){
        //gaussian kernel
        value = PARAM::SEARCH_RADIUS_H;
      }
      if(PARAM::KERNEL_FUNCTION==1){
        //cubic spline kernel
        value = 2.0;
      }
      return value * this->smth;
    }
    void setPos(const PS::F64vec& pos){
      this->pos = pos;
    }
  };
  class Grav{
  public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64    smth;
    PS::S64    id;
    PS::F64vec getPos() const{
      return this->pos;
    }
    PS::F64 getCharge(void) const{
      return this->mass;
    }
    void copyFromFP(const RealPtcl& rp){
      this->mass = rp.mass;
      this->smth = rp.smth;
      this->pos  = rp.pos;
      this->id   = rp.id;
    }
  };
}

#ifdef VARIABLE_SMOOTHING_LENGTH
#define Gather_or_Scatter Gather
#define Symmetry_or_Scatter Symmetry
#else
#define Gather_or_Scatter Scatter
#define Symmetry_or_Scatter Scatter
#endif
