/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <unistd.h>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 

  double interval = 0.05;
  int nsteps = 1.8 / interval;
  int nthermal = 1000;

  Input();

  for(int istep=0; istep < nsteps; istep++){

      for (int i=0; i<nspin; ++i)             //configurazione T=infinito
      {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      }
      
      for(int j=0; j<nthermal; j++){      //termalizzazione
        Move(metro);
      }

      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
      }
      ConfFinal(); //Write final configuration

      cout << "Measures at T = " << temp << " completed" << endl;

      temp -= interval;
      beta = 1.0/temp;
  }

  return 0;
}


void Input(void)
{


    ifstream ReadInput, Primes, input;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
  //Read input informations
    ReadInput.open("input.in");

    ReadInput >> restart;

    if(restart) input.open("seed.out");
    else input.open("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();

    ReadInput >> temp;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
      
    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;

    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;

    sleep(1);
    ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  ifstream ReadConf;

  if(restart){

    ReadConf.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
      ReadConf >> s[i];       //read configuration
    }
    ReadConf.close();

  }
  else{

    for (int i=0; i<nspin; ++i)             //configurazione T=infinito
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }

    // for(int i=0; i<nspin; i++){              //configurazione T=0
    //   s[i] = 1;
    // }
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;


    for (int i=0; i<nspin; ++i)             //configurazione T=infinito
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
    

}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis

      energy_old = Boltzmann(s[o], o);     //calcolo differenza di energia
      energy_new = Boltzmann(-s[o], o);

      double A = min(1., exp(- beta * (energy_new - energy_old)));   //acceptance probability

      double r = rnd.Rannyu();

      if(r < A){
        s[o] *= -1;
        accepted++;
      } 
    }
    else{ //Gibbs sampling

      energy_up = Boltzmann(1, o);
      energy_down = Boltzmann(-1, o);

      double prob = 1. / (1. + exp( -beta * (energy_up - energy_down)));   

      double r = rnd.Rannyu();

      if(r < prob){
        s[o] = -1;
      }
      else{
        s[o] = 1;
      } 

      accepted++;
    }
    attempted++;
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, c = 0.0, x = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);

     //c += u*u;

     x += s[i];

     m += s[i];
  }

  walker[iu] = u;     //divido dopo per nspin 
  walker[ic] = u*u;//beta*beta * (u*u/(double)nspin - pow(u/(double)nspin, 2));                 
  walker[ix] = beta * pow(x, 2) / (double)nspin;                                     
  walker[im] = m;               //divido dopo per nspin                    
                                                                           
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
  //ofstream Ene, Heat, Mag, Chi;
  const int wd=18;
   
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;

  //Ene.open("output_ene.dat",ios::app);
  //Heat.open("output_Cv.dat", ios::app);
  //Chi.open("output_chi.dat", ios::app);
  //Mag.open("output_mag.dat", ios::app);


  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);

  //stima_c = blk_av[ic]/blk_norm; //Heat capacity 
  stima_c = beta * beta * ( blk_av[ic]/blk_norm - blk_av[iu]/blk_norm * blk_av[iu]/blk_norm)/(double)nspin; 
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);

  stima_x = blk_av[ix]/blk_norm; //Susceptibility
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);

  stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);

  //Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  //Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  //Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  //Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;

  cout << "----------------------------" << endl << endl;

  if(iblk == nblk){

    ofstream Ene_curve, Heat_curve, Mag_curve, Chi_curve;

    Ene_curve.open("output_ene_curve.dat",ios::app);
    Heat_curve.open("output_Cv_curve.dat", ios::app);
    Chi_curve.open("output_chi_curve.dat", ios::app);
    Mag_curve.open("output_mag_curve.dat", ios::app);

    Ene_curve << setw(wd) << temp <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Heat_curve << setw(wd) << temp <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Chi_curve << setw(wd) << temp <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Mag_curve << setw(wd) << temp <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;

    cout << "Step completed" << endl;
    //sleep(1);

    Ene_curve.close();
    Heat_curve.close();
    Chi_curve.close();
    Mag_curve.close();

    //Ene.clear();
    //Heat.clear();
    //Chi.clear();
    //Mag.clear();
  }

  //Ene.close();
  //Heat.close();
  //Chi.close();
  //Mag.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
