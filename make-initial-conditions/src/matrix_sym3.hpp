#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>
#include"vector3.hpp"

namespace ParticleSimulator{
  //Matrix3 クラス 3x3の行列を記述
  //M = マトリクス, V=ベクトル, S=スカラー
  //M = M, M = S
  //M + M, M += M, +M, M - M, M -= M, -M
  //M * S, S * M, M *= S, M / S, M /= S
  //M * V (ベクトル積, Va = Mar * Vr rは和の規約で足し潰す)
  //M * M (MAij = MBir * MCjr rは和の規約に従って足し潰す)
  //上記の演算は定義・確認済み (2016/7/8 SugiuraK)
    template<class T>
    class Matrix3{
    public:
    //constructor
      T xx, xy, xz, yx, yy, yz, zx, zy, zz;
      Matrix3() : xx(T(0)), xy(T(0)), xz(T(0)), yx(T(0)), yy(T(0)), yz(T(0)), zx(T(0)), zy(T(0)), zz(T(0)) {} 
      Matrix3(const T _xx, const T _xy, const T _xz, 
	      const T _yx, const T _yy, const T _yz, 
	      const T _zx, const T _zy, const T _zz ) 
	: xx(_xx), xy(_xy), xz(_xz), yx(_yx), yy(_yy), yz(_yz), zx(_zx), zy(_zy), zz(_zz) {} 
      Matrix3(const T s) : xx(s), xy(s), xz(s), yx(s), yy(s), yz(s), zx(s), zy(s), zz(s) {}
      Matrix3(const Matrix3 & src) : xx(src.xx), xy(src.xy), xz(src.xz), yx(src.yx), yy(src.yy), yz(src.yz), zx(src.zx), zy(src.zy), zz(src.zz) {}

      const Matrix3 & operator = (const Matrix3 & rhs) {
	xx = rhs.xx;
	xy = rhs.xy;
	xz = rhs.xz;
	yx = rhs.yx;
	yy = rhs.yy;
	yz = rhs.yz;
	zx = rhs.zx;
	zy = rhs.zy;
	zz = rhs.zz;
	return (*this);
      }
      const Matrix3 & operator = (const T s) {
	xx = xy = xz = yx = yy = yz = zx = zy = zz = s;
	return (*this);
      }
      
      Matrix3 operator + (const Matrix3 & rhs) const {
	return Matrix3( xx + rhs.xx, xy + rhs.xy, xz + rhs.xz,
			yx + rhs.yx, yy + rhs.yy, yz + rhs.yz,
			zx + rhs.zx, zy + rhs.zy, zz + rhs.zz );
      }
      const Matrix3 & operator += (const Matrix3 & rhs) {
	(*this) = (*this) + rhs;
	return (*this);
      }
      const Matrix3 & operator + () const {
	return (* this);
      }
      Matrix3 operator - (const Matrix3 & rhs) const {
	return Matrix3( xx - rhs.xx, xy - rhs.xy, xz - rhs.xz,
			yx - rhs.yx, yy - rhs.yy, yz - rhs.yz,
			zx - rhs.zx, zy - rhs.zy, zz - rhs.zz );
      }
      const Matrix3 & operator -= (const Matrix3 & rhs) {
	(*this) = (*this) - rhs;
	return (*this);
      }
      const Matrix3 operator - () const {
	return Matrix3(-xx, -xy, -xz,
		       -yx, -yy, -yz,
		       -zx, -zy, -zz );
      }

      Matrix3 operator * (const T s) const{
	return Matrix3(xx * s, xy * s, xz * s,
		       yx * s, yy * s, yz * s,
		       zx * s, zy * s, zz * s );
      }
      friend Matrix3 operator * (const T s, const Matrix3 & mat){
	return (mat * s);
      }
      const Matrix3 & operator *= (const T s){
	(*this) = (*this) * s;
	return (*this);
      }
      Matrix3 operator / (const T s) const{
	return Matrix3(xx / s, xy / s, xz / s,
		       yx / s, yy / s, yz / s,
		       zx / s, zy / s, zz / s );
      }      
      const Matrix3 & operator /= (const T s){
	(*this) = (*this) / s;
	return (*this);
      }

      Vector3<T> operator * (const Vector3<T> v) const{
	return Vector3<T>(xx * v.x + xy * v.y + xz * v.z,
		       yx * v.x + yy * v.y + yz * v.z,
		       zx * v.x + zy * v.y + zz * v.z );
      }

      //MAij = MBir * MCjr rは和の規約に従って足し潰す
      Matrix3 operator * (const Matrix3 & rhs) const {
	return Matrix3( xx * rhs.xx + xy * rhs.xy + xz * rhs.xz, xx * rhs.yx + xy * rhs.yy + xz * rhs.yz, xx * rhs.zx + xy * rhs.zy + xz * rhs.zz,
			yx * rhs.xx + yy * rhs.xy + yz * rhs.xz, yx * rhs.yx + yy * rhs.yy + yz * rhs.yz, yx * rhs.zx + yy * rhs.zy + yz * rhs.zz,
			zx * rhs.xx + zy * rhs.xy + zz * rhs.xz, zx * rhs.yx + zy * rhs.yy + zz * rhs.yz, zx * rhs.zx + zy * rhs.zy + zz * rhs.zz );
      }

      T getTrace() const {
	return (xx + yy + zz);
      }

      T getJ2() const {
	return( 0.5 *(xx*xx + xy*xy + xz*xz + yx*yx + yy*yy + yz*yz + zx*zx + zy*zy + zz*zz) );
      }

      template <typename U>
      operator Matrix3<U> () const {
	return Matrix3<U>( static_cast<U>(xx), static_cast<U>(xy), static_cast<U>(xz),
			   static_cast<U>(yx), static_cast<U>(yy), static_cast<U>(yz),
			   static_cast<U>(zx), static_cast<U>(zy), static_cast<U>(zz) );
      }
      
      friend std::ostream & operator << (std::ostream & c, const Matrix3 & mat){
	c<<mat.xx<<"   "<<mat.xy<<"    "<<mat.xz<<std::endl;
	c<<mat.xy<<"   "<<mat.yy<<"    "<<mat.yz<<std::endl;
	c<<mat.xz<<"   "<<mat.yz<<"    "<<mat.zz<<std::endl;
	return c;
      }
    };

    template<class T>
    class MatrixSym3{
    public:
        //constructor
        T xx, yy, zz, xy, xz, yz;
        MatrixSym3() : xx(T(0)), yy(T(0)), zz(T(0)), xy(T(0)), xz(T(0)), yz(T(0)) {} 
        MatrixSym3(const T _xx, const T _yy, const T _zz, 
                   const T _xy, const T _xz, const T _yz ) 
            : xx(_xx), yy(_yy), zz(_zz), xy(_xy), xz(_xz), yz(_yz) {} 
        MatrixSym3(const T s) : xx(s), yy(s), zz(s), xy(s), xz(s), yz(s) {}
        MatrixSym3(const MatrixSym3 & src) : xx(src.xx), yy(src.yy), zz(src.zz), 
                                             xy(src.xy), xz(src.xz), yz(src.yz) {}

        const MatrixSym3 & operator = (const MatrixSym3 & rhs) {
            xx = rhs.xx;
            yy = rhs.yy;
            zz = rhs.zz;
            xy = rhs.xy;
            xz = rhs.xz;
            yz = rhs.yz;
            return (*this);
        }
        const MatrixSym3 & operator = (const T s) {
            xx = yy = zz = xy = xz = yz = s;
            return (*this);
        }

        MatrixSym3 operator + (const MatrixSym3 & rhs) const {
            return MatrixSym3( xx + rhs.xx, yy + rhs.yy, zz + rhs.zz,
                               xy + rhs.xy, xz + rhs.xz, yz + rhs.yz);
        }
        const MatrixSym3 & operator += (const MatrixSym3 & rhs) const {
            (*this) = (*this) + rhs;
            return (*this);
        }
        MatrixSym3 operator - (const MatrixSym3 & rhs) const {
            return MatrixSym3( xx - rhs.xx, yy - rhs.yy, zz - rhs.zz,
                               xy - rhs.xy, xz - rhs.xz, yz - rhs.yz);
        }
        const MatrixSym3 & operator -= (const MatrixSym3 & rhs) const {
            (*this) = (*this) - rhs;
            return (*this);
        }

        T getTrace() const {
            return (xx + yy + zz);
        }

        template <typename U>
        operator MatrixSym3<U> () const {
            return MatrixSym3<U>( static_cast<U>(xx), static_cast<U>(yy), static_cast<U>(zz),
                                  static_cast<U>(xy), static_cast<U>(xz), static_cast<U>(yz) );
        }

        friend std::ostream & operator << (std::ostream & c, const MatrixSym3 & mat){
            c<<mat.xx<<"   "<<mat.xy<<"    "<<mat.xz<<std::endl;
            c<<mat.xy<<"   "<<mat.yy<<"    "<<mat.yz<<std::endl;
            c<<mat.xz<<"   "<<mat.yz<<"    "<<mat.zz<<std::endl;
            return c;
        }
    };
}
