#ifndef __RWclasses__
#define __RWclasses__

#include "random.h"
#include "eserc2functions.h"

using namespace std;

class path{

    public:

        path(Random& r, int n = 1E2) : paths(n, vector<double>(3,0)), num{n}, norm2(n,0), rnd{r} {}

        void StepDiscreto(){       //Fa avanzare tutti i walker di uno step discreto

            int r = 0;

            for(int i=0; i < num; i++){

                r = floor(rnd.Rannyu(0.,6.));  //0,1,2 significa avanti in x,y,z rispettivamente. 3,4,5 indietro.

                if(r<3){
                    paths[i][r]++;
                }
                else{
                    paths[i][r-3]--;
                }

                norm2[i] = pow(paths[i][0],2) + pow(paths[i][1],2) + pow(paths[i][2],2);    //calcolo norma
            }
        }

        void StepContinuo(){         //Fa avanzare tutti i walker di uno step continuo

            double r = 0.;
            double theta = 0.;
            double phi = 0.;

            for(int i=0; i < num; i++){

                r = rnd.Rannyu();

                theta = acos(1. - 2.*r);
                phi = rnd.Rannyu(0., 2*M_PI);

                paths[i][0] += sin(theta) * cos(phi);
                paths[i][1] += sin(theta) * sin(phi);
                paths[i][2] += cos(theta);

                norm2[i] = pow(paths[i][0],2) + pow(paths[i][1],2) + pow(paths[i][2],2);    //calcolo norma
            }
        }

        vector<double> GetNorm2(){
            return norm2;
        }

        int GetNum(){
            return num;
        }

    public:
        vector<vector<double>> paths;
        int num;
        vector<double> norm2;
        Random rnd;
};

class BlockAverageRW{

    public:

        BlockAverageRW(path& fun, int n = 2, int b = 20):
            n_props{n}, nblk{b}, blk_av(n_props), glob_av(n_props), blk_av2(n_props), glob_av2(n_props), walker(n_props), p{fun} {
                blk_norm = 0.;
                stimaP = 0.;
                stimaEr = 0.;
                errore = 0.;
            }
        
        //////////////////  Metodi  ////////////////////////

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

        void Measure(int i){

            walker[ipath] = p.norm2[i];
            walker[ierr] = sqrt(p.norm2[i]);
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk, ofstream& output){ //Print results for current block
                   
            const int wd=12;

            stimaP = blk_av[ipath]/blk_norm;     //calcolo media di r^2
            glob_av[ipath] += stimaP;
            glob_av2[ipath] += stimaP*stimaP;

            stimaEr = blk_av[ierr]/blk_norm;     //calcolo media di sqrt(r^2)  per il calcolo dell'errore
            glob_av[ierr] += stimaEr;
            glob_av2[ierr] += stimaEr*stimaEr;

            if(iblk == nblk){      //stampo solo l'ultimo valore, quello del blocco finale.
                
                //errore = error(glob_av[ipath],glob_av2[ipath],iblk);  //errore su r^2
                errore = error(glob_av[ierr],glob_av2[ierr],iblk);       //errore su sqrt(r^2)

                output << setw(wd) << iblk <<  setw(wd) << sqrt(stimaP) << setw(wd) << sqrt(glob_av[ipath]/(double)iblk) << setw(wd) << errore << endl;
                //output << setw(wd) << iblk <<  setw(wd) << sqrt(stimaP) << setw(wd) << sqrt(glob_av[ipath]/(double)iblk) << setw(wd) << errore / (2 * sqrt(glob_av[ipath]/(double)iblk)) << endl;
                // con propagazione errori, viene piu' piccolo di un fattore 10 rispetto a errore su r^2. errore su sqrt(r^2) e' dello stesso ordine di grandezza.
            }

        }


    public:
    
        int n_props, nblk;
        const int ipath = 0;
        const int ierr = 1;

        vector<double>  blk_av, glob_av, blk_av2, glob_av2;
        vector<double>  walker;
        double stimaP, stimaEr;
        double errore;
        double blk_norm;

        path p;

};


#endif 