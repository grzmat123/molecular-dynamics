#ifndef _dynamika_h_
#define _dynamika_h_
#include <iostream>
#include <cmath>
#include <conio.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>

using namespace std;

const double k = 0.00831;

class Krysztal{
      
      public:
      
      double b0[3];
      double b1[3];
      double b2[3];       
      double **polozenia;
      double **energie_kin;
      double **pedy;
      double **sily_scianki;
      double **sily_pot;
      double potencjal;
      double cisnienie;
      double T_aver;
      double P_aver;
      double H_aver;
      int n;
      int N;
      double a;
      double To;
      double m;
      double L;
      double f;
      double R;
      double e;
      double H;
      double T;
      Krysztal(double nn = 0, double aa = 0, double Too = 0, double mm = 0, double ll = 0, double ff = 0, double rr = 0, double ee = 0);
      ~Krysztal(){}
      friend void generuj_krysztal(Krysztal arg);
      friend void zapisz_sily(Krysztal arg);
      void normalizuj();
      void eliminuj_ruch_srodka();
      void licz_sily();
      void licz_hamiltonian();
      void licz_temperature();
      void licz_cisnienie();
      void symulacja(double tau, double So, double Sd, int Sout, int Sxyz);
      };
      
double *wczytaj_dane(double *t);
void generuj_krysztal(Krysztal arg);
void zapisz_sily(Krysztal arg);
void ustal_a();
#endif
