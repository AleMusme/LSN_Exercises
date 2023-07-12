#ifndef __Buffonclasses__
#define __Buffonclasses__

#include "random.h"
#include "eserc1functions.h"

using namespace std;

class pi{

    public:

        pi(){}

        pi(Random& r, double l = 0.75, double d = 1.) : rnd{r}{
            accepted = 0;
            attempted = 0;
            length = l;
            distance = d;
        }

        double Eval(void){
            return 2. * length * (double) attempted / ((double) accepted * distance);    //PROBLEMA: se metto accepted=attempted=0 rischio infinito al primo passo.
                                                                                         //Va bene mettere =1 inizialmente entrambi?
        }


        void Step(){    //Le linee orizzontali hanno distanza d tra loro. 
                        //Si estrae un numero tra 0 e d che è un estremo dell'ago. si estrae un angolo uniformemente tra 0 e 2pi (vedere random.cpp)
                        //e si determina il secondo estremo. se il secondo punto ha ordinata maggiore di d, o minore di 0, allora interseca una linea.
        
            double estremo1 = rnd.Rannyu(0., distance);
            double estremo2 = estremo1 + length * sin(rnd.UnifTheta());

            if(attempted==0){
                accepted++;
            }
            else{
                if(estremo2 > distance or estremo2 < 0){ 
                    accepted++;
                }
            }
            attempted++;
        }
        
    public:
        int accepted;
        int attempted;
        double length;
        double distance;

        Random rnd;

};

class BlockAveragePi{

    public:

        BlockAveragePi(pi& fun, int b = 20, int n = 1) :
            n_props{n}, nblk{b}, blk_av(n_props), glob_av(n_props), glob_av2(n_props), walker(n_props), function{fun} {
                blk_norm = 0.;
                stima = 0.;
                errore = 0.;
            }

        void Reset(int iblk){ //Reset block averages
    
            if(iblk == 1){

                for(int i=0; i<n_props; ++i)
                {
                    glob_av[i] = 0;
                    glob_av2[i] = 0;
                }
            }

            for(int i=0; i<n_props; ++i){

                blk_av[i] = 0;
            }

            blk_norm = 0;
        }

        void Measure(void){

            function.Step();
            walker[iav] = function.Eval();
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block
        
            ofstream Pi; 
            const int wd=18;
            
            //cout << "Block number " << iblk << endl;

            Pi.open("output_pi.dat",ios::app);

            stima = blk_av[iav]/blk_norm; 
            glob_av[iav] += stima;
            glob_av2[iav] += stima*stima;
            errore = error(glob_av[iav],glob_av2[iav],iblk);

            Pi << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iav]/(double)iblk << setw(wd) << errore << endl;

            Pi.close();
        }


    public:
    
        int n_props, nblk;
        const int iav=0;
        const int ierr=1;

        vector<double>  blk_av, glob_av, glob_av2;
        vector<double> walker;
        double stima;
        double errore;
        double blk_norm;

        pi function;

};


#endif 