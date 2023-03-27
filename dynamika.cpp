#include "dynamika.h"


Krysztal::Krysztal(double nn, double aa, double Too, double mm, double ll, double ff, double rr, double ee){
             n = nn;
             N = n*n*n;
             a = aa;
             To = Too;
             m = mm;
             L = ll;
             f = ff;
             R = rr;
             e = ee;
             H = 0;
             T = To;
             T_aver = 0;
             P_aver = 0;
             H_aver = 0;
             b0[0] = a;
             b0[1] = 0;
             b0[2] = 0;
             b1[0] = a/2;
             b1[1] = sqrt(3)*a/2;
             b1[2] = 0;
             b2[0] = a/2;
             b2[1] = sqrt(3)*a/6;
             b2[2] = a*sqrt(2)/sqrt(3);
             polozenia = new double *[N];
             energie_kin = new double *[N];
             pedy = new double *[N];
             sily_scianki = new double *[N];
             sily_pot = new double *[N];
              for (int i = 0; i < N ; i++)
              {
               polozenia[i] = new double[3];
               energie_kin[i] = new double[3];
               pedy[i] = new double[3]; 
               sily_scianki[i] = new double[3];
               sily_pot[i] = new double[3];
               }
               for (int i0 = 0 ; i0 < n ; i0 ++)
               {
                for (int i1 = 0 ; i1 < n ; i1 ++)
                {
                 for (int i2 = 0 ; i2 < n ; i2 ++)
                 {
                 int i;
                 i = i0 + i1*n + i2*n*n;
        
                 polozenia[i][0] = (i0 - (n-1)/2)*b0[0] + (i1 - (n-1)/2)*b1[0] + (i2 - (n-1)/2)*b2[0];
                 polozenia[i][1] = (i1 - (n-1)/2)*b1[1] + (i2 - (n-1)/2)*b2[1];
                 polozenia[i][2] = (i2 - (n-1)/2)*b2[2];  
                
                double loteriai0 =( rand() % 100 ) + 1;
                       double lambdai0 = loteriai0/100;
                              double lni0 = log( lambdai0 );
                double loteriai1 =( rand() % 100 ) + 1;
                       double lambdai1 = loteriai1/100;
                              double lni1 = log( lambdai1 );
                double loteriai2 =( rand() % 100 ) + 1;
                       double lambdai2 = loteriai2/100;
                              double lni2 = log( lambdai2 );
                energie_kin[i][0] = -0.5 * k * To * lni0;
                energie_kin[i][1] = -0.5 * k * To * lni1;
                energie_kin[i][2] = -0.5 * k * To * lni2;
       }   
      } 
      }
      normalizuj();
          for ( int i = 0 ; i < N ; i++ ) // petla do obliczania poczatkowych wartosci pedu
    {
        double loteriapi0 =( rand() % 2 ) + 1;
               double znakpi0 = (loteriapi0 == 2) ? 1 : -1;
        double loteriapi1 =( rand() % 2 ) + 1;
               double znakpi1 = (loteriapi1 == 2) ? 1 : -1;
        double loteriapi2 =( rand() % 2 ) + 1;
               double znakpi2 = (loteriapi2 == 2) ? 1 : -1;
        pedy[i][0] = znakpi0 * sqrt(2*m*energie_kin[i][0]); 
        pedy[i][1] = znakpi1 * sqrt(2*m*energie_kin[i][1]);
        pedy[i][2] = znakpi2 * sqrt(2*m*energie_kin[i][2]);
    }
    eliminuj_ruch_srodka(); 
    licz_sily();
    licz_cisnienie();
    licz_hamiltonian();                    
    }

void Krysztal::normalizuj(){
               double sum = 0;
               if(To != 0){
               for ( int i = 0 ; i < N ; i++ )
               {
               sum = sum + energie_kin[i][0] + energie_kin[i][1] + energie_kin[i][2]; 
               }
               double srednia = sum/(3*N);
               double przelicznik = (0.5*k*To)/srednia; // obliczenie wsp. do normalizacji
               for ( int i = 0 ; i < N ; i++ ) // petla do obliczania znorm. en. kin.
               {
               energie_kin[i][0] = przelicznik * energie_kin[i][0]; 
               energie_kin[i][1] = przelicznik * energie_kin[i][1];
               energie_kin[i][2] = przelicznik * energie_kin[i][2];
               }}
               else
               {
               for ( int i = 0 ; i < N ; i++ ) // petla do obliczania znorm. en. kin.
               {
               energie_kin[i][0] = 0; 
               energie_kin[i][1] = 0;
               energie_kin[i][2] = 0;
               } 
                   }
               }
  
void Krysztal::eliminuj_ruch_srodka(){
     double sump0 = 0, sump1 = 0, sump2 = 0;
     for ( int i = 0 ; i < N ; i++ ) // petla pomocnicza do eliminacji ruchu sr. masy
     {
        sump0 = sump0 + pedy[i][0];
        sump1 = sump1 + pedy[i][1];
        sump2 = sump2 + pedy[i][2];
        }
        for ( int i = 0 ; i < N ; i++ ) // petla do obliczania pedow przy wyeliminowanym ruchu sr. masy
        {
        pedy[i][0] = pedy[i][0] - sump0/N; 
        pedy[i][1] = pedy[i][1] - sump1/N;
        pedy[i][2] = pedy[i][2] - sump2/N;
        }
     } 

void Krysztal::licz_sily(){
         potencjal = 0;
         for (int i = 0; i < N; i++){
        
        sily_pot[i][0] = 0;
        sily_pot[i][1] = 0;
        sily_pot[i][2] = 0;
        
        }
     for ( int i = 0 ; i < N  ; i++ ) //petla do akumulacji potencjalow i sil odpychania od scianek oraz cisnienia chwilowego
    {
        double modri = 0;
         modri = sqrt(polozenia[i][0]*polozenia[i][0] + polozenia[i][1]*polozenia[i][1] + polozenia[i][2]*polozenia[i][2]);
        double VSri = 0;
        if (modri < L)
        {
        VSri = 0;          
        }
        else
        {
        VSri = 0.5 * f * (modri - L) * (modri - L);    
        }

        potencjal = potencjal + VSri;

        if (modri < L)
        {
        sily_scianki[i][0] =  0;
        sily_scianki[i][1] =  0;
        sily_scianki[i][2] =  0;      
        }
        else
        {
        sily_scianki[i][0] =  f * (L - modri) * polozenia[i][0]/modri;
        sily_scianki[i][1] =  f * (L - modri) * polozenia[i][1]/modri; 
        sily_scianki[i][2] =  f * (L - modri) * polozenia[i][2]/modri;     
        }

                      if ( i > 0 )
                      {
                         for ( int j = 0 ; j < i ; j ++)
                        {
                           double modrij = sqrt((polozenia[i][0]-polozenia[j][0])*(polozenia[i][0]-polozenia[j][0])
                           + (polozenia[i][1]-polozenia[j][1])*(polozenia[i][1]-polozenia[j][1])
                            + (polozenia[i][2]-polozenia[j][2])*(polozenia[i][2]-polozenia[j][2]));

                           double Vprij = e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6));
                           
                           potencjal = potencjal + Vprij;
                           sily_pot[i][0] = sily_pot[i][0] + 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][0] - polozenia[j][0])/(modrij * modrij);
                           sily_pot[i][1] = sily_pot[i][1] + 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][1] - polozenia[j][1])/(modrij * modrij);
                           sily_pot[i][2] = sily_pot[i][2] + 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][2] - polozenia[j][2])/(modrij * modrij);
                           sily_pot[j][0] = sily_pot[j][0] - 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][0] - polozenia[j][0])/(modrij * modrij);
                           sily_pot[j][1] = sily_pot[j][1] - 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][1] - polozenia[j][1])/(modrij * modrij);
                           sily_pot[j][2] = sily_pot[j][2] - 12 * e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6)) * (polozenia[i][2] - polozenia[j][2])/(modrij * modrij);                             
                         }                    
                      }
    }
     }

void Krysztal::licz_hamiltonian(){
     H = 0;
     for (int i = 0; i < N ; i++){
         H = H + (pedy[i][0]*pedy[i][0] + pedy[i][1]*pedy[i][1] + pedy[i][2]*pedy[i][2])*0.5/m;
         }
     H = H + potencjal;
     }

void Krysztal::licz_temperature(){
     licz_hamiltonian();
     double E = H - potencjal;
     T = 2*E/(3*N*k);
     }

void Krysztal::licz_cisnienie(){
          cisnienie = 0;                      
          for (int i = 0; i < N ; i++){
          double modFiS = sqrt(sily_scianki[i][0]*sily_scianki[i][0] 
          + sily_scianki[i][1]*sily_scianki[i][1] 
          + sily_scianki[i][2]*sily_scianki[i][2]);
          cisnienie = cisnienie + modFiS/(4.*M_PI*L*L);
                         } 
     }
  
void Krysztal::symulacja(double tau, double So, double Sd, int Sout, int Sxyz){
     fstream output;
     output.open ("wartosci_chwilowe.txt", ios::out);
     if( output.good() == true )
    {
    cout << "Zapisuje wartosci chwilowe H, V, P, T do pliku wartosci_chwilowe!" << endl;
    output<<"t H V T P"<<endl;
    }
    fstream aver;
     aver.open ("symulacja.xyz", ios::out);
     if( aver.good() == true )
    {
    cout << "Zapisuje klatki symulacji do pliku symulacja!" << endl;
    }
    else{cout<<"Brak dostepu do symulacja!"<<endl;}
     for( int s = 1; s < So + Sd ; s++){
        for (int i = 0; i < N; i++){
            pedy[i][0] = pedy[i][0] + 0.5*(sily_pot[i][0]+sily_scianki[i][0])*tau;
            pedy[i][1] = pedy[i][1] + 0.5*(sily_pot[i][1]+sily_scianki[i][1])*tau;
            pedy[i][2] = pedy[i][2] + 0.5*(sily_pot[i][2]+sily_scianki[i][2])*tau;
            polozenia[i][0] = polozenia[i][0] + pedy[i][0]*tau/m;
            polozenia[i][1] = polozenia[i][1] + pedy[i][1]*tau/m;
            polozenia[i][2] = polozenia[i][2] + pedy[i][2]*tau/m;
            }
            licz_sily();
            for (int i = 0; i < N; i++){
            pedy[i][0] = pedy[i][0] + 0.5*(sily_pot[i][0]+sily_scianki[i][0])*tau;
            pedy[i][1] = pedy[i][1] + 0.5*(sily_pot[i][1]+sily_scianki[i][1])*tau;
            pedy[i][2] = pedy[i][2] + 0.5*(sily_pot[i][2]+sily_scianki[i][2])*tau;
            }
            licz_temperature();
            licz_hamiltonian();
            licz_cisnienie();      
        if (s%Sout == 0) {output<<s*tau<<" "<<H<<" "<<potencjal<<" "<<T<<" "<<cisnienie<<endl;};
        if (s%Sxyz == 0){aver<<N<<endl<<"#"<<endl;
                   for (int i = 0; i < N; i++){
                   aver<<"Ar"<<" "<<polozenia[i][0]<<" "<<polozenia[i][1]<<" "<<polozenia[i][2]<<endl;
                   
                   }                  
                   };     
        if (s >= So) {T_aver = T_aver + T;
                     P_aver = P_aver + cisnienie;
                     H_aver = H_aver + H;}            
        }
        T_aver = T_aver/Sd;
        P_aver = P_aver/Sd;
        H_aver = H_aver/Sd;
        output<<"H_aver = "<<H_aver<<" P_aver = "<<P_aver<<" T_aver = "<<T_aver<<endl;
        output.close();
        aver.close();
     }
               
double *wczytaj_dane(double *t){
     char cos[200];
     fstream param ("param.txt", ios::in);
         for(int i = 0 ; i < 13 ; i ++){
            param >> t[i];
            param >> cos;
            };
     return t;    
     }
     
void generuj_krysztal(Krysztal arg){
         fstream argon;
    argon.open ("argon.xyz" , ios::out);
    if( argon.good() == true )
    {
    cout << "Zapisuje poczatkowa postac krysztalu do pliku argon!" << endl;
    argon<<arg.N<<endl<<endl;
    for (int i = 0 ; i < arg.N ; i++){
    argon<<"Ar "<<arg.polozenia[i][0]<<" "<<arg.polozenia[i][1]<<" "<<arg.polozenia[i][2]<<endl;}
    } else cout << "Dostep do pliku argon zostal zabroniony!" << endl;
    argon.close();
     }
void zapisz_sily(Krysztal arg){
              fstream sily;
    sily.open ("sily.txt" , ios::out);
    if( sily.good() == true )
    {
    cout << "Zapisuje sily do pliku sily!" << endl;
    sily<<arg.N<<endl<<endl;
    for (int i = 0 ; i < arg.N ; i++){
    sily<<"Ar "<<arg.sily_scianki[i][0] + arg.sily_pot[i][0]<<" "<<arg.sily_scianki[i][1] + arg.sily_pot[i][1]<<" "<<arg.sily_scianki[i][2] + arg.sily_pot[i][2]<<endl;}
    } else cout << "Dostep do pliku sily zostal zabroniony!" << endl;
    sily<<arg.potencjal<<endl;
    sily<<arg.cisnienie<<endl;
    sily.close();
    }

void ustal_a(){
     double a;
     int n = 5;
     double b0[3];
     double b1[3];
     double b2[3];
     double V;
     double e = 1;
     double R = 0.38;
     double **polozenia;
     fstream energia;
     energia.open ("energia.txt", ios::out);
     energia<<"a "<<"V"<<endl;
     for (int i = 0; i < 100; i++){
         V = 0;
         delete polozenia;
         a = 0.3 + 0.001*i;
             b0[0] = a;
             b0[1] = 0;
             b0[2] = 0;
             b1[0] = a/2;
             b1[1] = sqrt(3)*a/2;
             b1[2] = 0;
             b2[0] = a/2;
             b2[1] = sqrt(3)*a/6;
             b2[2] = a*sqrt(2)/sqrt(3);
             polozenia = new double *[125];
              for (int i = 0; i < 125; i ++)
              {
               polozenia[i] = new double[3];
               }
               for (int i0 = 0 ; i0 < 5; i0 ++)
               {
                for (int i1 = 0 ; i1 < 5; i1 ++)
                {
                 for (int i2 = 0 ; i2 < 5; i2 ++)
                 {
                 int i;
                 i = i0 + i1*5+ i2*5*5;
                 polozenia[i][0] = (i0 - (n-1)/2)*b0[0] + (i1 - (n-1)/2)*b1[0] + (i2 - (n-1)/2)*b2[0];
                 polozenia[i][1] = (i1 - (n-1)/2)*b1[1] + (i2 - (n-1)/2)*b2[1];
                 polozenia[i][2] = (i2 - (n-1)/2)*b2[2];
                 }
                 }
                 }
                 for (int i = 0; i < 125; i++){
                 if ( i > 0 )
                      {
                         for ( int j = 0 ; j < i ; j ++)
                        {
                           double modrij = sqrt((polozenia[i][0]-polozenia[j][0])*(polozenia[i][0]-polozenia[j][0])
                           + (polozenia[i][1]-polozenia[j][1])*(polozenia[i][1]-polozenia[j][1])
                            + (polozenia[i][2]-polozenia[j][2])*(polozenia[i][2]-polozenia[j][2]));

                           double Vprij = e * (pow(R/modrij, 12) - 2*pow(R/modrij, 6));
                           
                           V = V + Vprij;
                           }
                           }
                           }
        energia<<a<<" "<<V<<endl;
         
         }
  
      delete  polozenia;
      energia.close();
     }
