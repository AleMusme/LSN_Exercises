#ifndef __eserc1classes__
#define __eserc1classes__

#include "random.h"
#include "eserc2functions.h"

using namespace std;

class f{

    public:

        double Eval(double x){
            return  M_PI/2. * cos(M_PI * x / 2.);
        }

        double EvalDistro(double x){
            return M_PI/4. * cos(M_PI * x / 2.) / (1. - x);
        }

};

class BlockAverage{

    public:

        BlockAverage(f& fun, int n = 2, int b = 20):
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

        void Measure(double x, double y){

            walker[inor] = function.Eval(x);
            walker[idis] = function.EvalDistro(y);
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block
        
            ofstream Normal, Distro; 
            const int wd=18;
            
            cout << "Block number " << iblk << endl;

            Normal.open("output_normal.dat",ios::app);
            Distro.open("output_distro.dat", ios::app);

            stima = blk_av[inor]/blk_norm; // normal integral
            glob_av[inor] += stima;
            glob_av2[inor] += stima*stima;
            errore = error(glob_av[inor],glob_av2[inor],iblk);

            Normal << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[inor]/(double)iblk << setw(wd) << errore << endl;

            stima = blk_av[idis]/blk_norm; // distro integral
            glob_av[idis] += stima;
            glob_av2[idis] += stima*stima;
            errore = error(glob_av[idis],glob_av2[idis],iblk);

            Distro << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[idis]/(double)iblk << setw(wd) << errore << endl;

            cout << "----------------------------" << endl << endl;

            Normal.close();
            Distro.close();

        }


    public:
    
        int n_props, nblk;
        const int inor=0;
        const int idis=1;

        vector<double>  blk_av, glob_av, glob_av2;
        vector<double> walker;
        double stima;
        double errore;
        double blk_norm;

        f function;

};


#endif 