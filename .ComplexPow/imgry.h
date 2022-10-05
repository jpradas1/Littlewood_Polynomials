#ifndef _IMGRY_H_
#define _IMGRY_H_

#include <iostream>
#include <stdio.h>
#include <math.h>

template <class T>
class Imaggi
{
    public:
        Imaggi ();
        Imaggi (const T& x, const T& y);
        Imaggi (const Imaggi<T>& z);

        ~Imaggi ();

        T GetRealPart (); T GetImagPart ();
        T GetPhase (const Imaggi<T> z);
        T GetAmplitud (const Imaggi<T> z);


        //operation on imaginary unit

        bool operator== (const Imaggi<T>& z);

        template <class R> friend Imaggi<R> operator+ (const Imaggi<R>& cxl, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator+ (const R& l, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator+ (const Imaggi<R>& cxl, const R& r);

        template <class R> friend Imaggi<R> operator- (const Imaggi<R>& cxl,const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator- (const R& l,const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator- (const Imaggi<R>& cxl, const R& r);

        template <class R> friend Imaggi<R> operator* (const Imaggi<R>& cxl, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator* (const R& l, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator* (const Imaggi<R>& cxl, const R& r);

        template <class R> friend Imaggi<R> operator/ (const Imaggi<R>& cxl, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator/ (const R& l, const Imaggi<R>& cxr);
        template <class R> friend Imaggi<R> operator/ (const Imaggi<R>& cxr, const R& r);

        template <class R> friend Imaggi<R> operator^ (const Imaggi<R>& z, const int n);

        template <class R> friend Imaggi<R> Conjugate (const Imaggi<R>& z);

    private:
        T *img;
        T real_p, imag_p;

};

// #######################   Constructors   #######################

template <class T>
Imaggi<T>::Imaggi (){
    real_p = 0.0;
    imag_p = 1.0;
    img = new T [2];
    img[0] = real_p; img[1] = imag_p;
}

template <class T>
Imaggi<T>::Imaggi (const T& x, const T& y){
    real_p = x;
    imag_p = y;
    img = new T [2];
    img[0] = real_p; img[1] = imag_p;
}

template <class T>
Imaggi<T>::Imaggi (const Imaggi<T>& z){
    real_p = z.real_p;
    imag_p = z.imag_p;
    img = new T [2];
    img[0] = real_p; img[1] = imag_p;
}

template <class T>
Imaggi<T>::~Imaggi (){
    if(img != nullptr)
        delete [] img;
}

// ####################### Some elemental functions #######################

template <class T>
T Imaggi<T>::GetRealPart(){

    return real_p;
}

template <class T>
T Imaggi<T>::GetImagPart(){

    return imag_p;
}

template <class T>
T Imaggi<T>::GetAmplitud(const Imaggi<T> z){
    T x = z.real_p; T y = z.imag_p;
    T r2 = x*x+y*y;
    T r = std::sqrt(r2);
    return r;
}

template <class T>
T Imaggi<T>::GetPhase(const Imaggi<T> z){
    T x = z.real_p; T y = z.imag_p;
    T phi = atan (x/y); //Radians
    return phi;
}

// ####################### Operations on Imaginaries #######################

// SUM
template <class T>
Imaggi<T> operator+ (const Imaggi<T>& cxl, const Imaggi<T>& cxr){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T xr = cxr.real_p;  T yr = cxr.imag_p;
    T X = xl + xr; T Y = yl + yr;

    Imaggi<T> result(X, Y);
    return result;
}

template <class T>
Imaggi<T> operator+ (const T& l, const Imaggi<T>& cxr){
    T xr = cxr.real_p;  T yr = cxr.imag_p;
    T X = l + xr;

    Imaggi<T> result(X, yr);
    return result;
}

template <class T>
Imaggi<T> operator+ (const Imaggi<T>& cxl, const T& r){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T X = xl + r;

    Imaggi<T> result(X, yl);
    return result;
}

//SUBTRACTION
template <class T>
Imaggi<T> operator- (const Imaggi<T>& cxl, const Imaggi<T>& cxr){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T xr = cxr.real_p;  T yr = cxr.imag_p;
    T X = xl - xr; T Y = yl - yr;

    Imaggi<T> result(X, Y);
    return result;
}

template <class T>
Imaggi<T> operator- (const T& l, const Imaggi<T>& cxr){
    T xr = cxr.real_p;  T yr = cxr.imag_p;
    T X = l - xr;

    Imaggi<T> result(X, yr);
    return result;
}

template <class T>
Imaggi<T> operator- (const Imaggi<T>& cxl, const T& r){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T X = xl - r;

    Imaggi<T> result(X, yl);
    return result;
}

//MULTIPLICATION
template <class T>
Imaggi<T> operator* (const Imaggi<T>& cxl, const Imaggi<T>& cxr){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T xr = cxr.real_p;  T yr = cxr.imag_p;

    T X = xl*xr - yl*yr;
    T Y = xl*yr + yl*xr;
    Imaggi<T> result(X, Y);
    return result;
}

template <class T>
Imaggi<T> operator* (const T& l, const Imaggi<T>& cxr){
    T xr = cxr.real_p;  T yr = cxr.imag_p;
    T X = l * xr; T Y = l * yr;

    Imaggi<T> result(X, Y);
    return result;
}

template <class T>
Imaggi<T> operator* (const Imaggi<T>& cxl, const T& r){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T X = xl * r; T Y = yl * r;

    Imaggi<T> result(X, Y);
    return result;
}

//DIVISION
template <class T>
Imaggi<T> operator/ (const Imaggi<T>& cxl, const Imaggi<T>& cxr){
    // cxl/cxr
    T xr = cxr.real_p; T yr = cxr.imag_p;
    T R2 = xr*xr + yr*yr;
    Imaggi<T> conj = Conjugate(cxr);
    Imaggi<T> result = (cxl * conj)/(R2);
    return result;
}

template <class T>
Imaggi<T> operator/ (const T& l, const Imaggi<T>& cxr){
    // l/cxr
    T xr = cxr.real_p; T yr = cxr.imag_p;
    T R2 = xr*xr + yr*yr;
    Imaggi<T> conj = Conjugate(cxr);
    Imaggi<T> result = (l * conj)/(R2);
    return result;
}

template <class T>
Imaggi<T> operator/ (const Imaggi<T>& cxl, const T& r){
    T xl = cxl.real_p;  T yl = cxl.imag_p;
    T X = xl / r; T Y = yl / r;

    Imaggi<T> result(X, Y);
    return result;
}

//EXPONENTIATION
template <class T>
Imaggi<T> operator^ (const Imaggi<T>& z, const int n){
    Imaggi<T> partial_res = z;
    for(int ii = 0; ii < (n-1); ii++){
        partial_res  = partial_res * z;
    }
    Imaggi<T> result(partial_res);
    return result;
}

//CONJUGATE
template <class T>
Imaggi<T> Conjugate (const Imaggi<T>& z){
    T real = z.real_p; T imag = z.imag_p;
    T I = -imag;

    Imaggi<T> result(real, I);
    return result;
}

//EQUALITY
template <class T>
bool Imaggi<T>::operator== (const Imaggi<T>& z){

    if((this-> real_p != z.real_p) && (this-> imag_p != z.imag_p))
        return false;
}

#endif
