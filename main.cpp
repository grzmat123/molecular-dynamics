#include "dynamika.h"

int main(){
double *temp = new double[13];
 srand( time( NULL ) );
 wczytaj_dane(temp);
 Krysztal Argon1(temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7]);
 generuj_krysztal(Argon1);
 Argon1.normalizuj();
 //zapisz_sily(Argon1);
 Argon1.symulacja(temp[8], temp[9], temp[10], temp[11], temp[12]);
 delete temp;
 //getch();
 return 0;   
}
