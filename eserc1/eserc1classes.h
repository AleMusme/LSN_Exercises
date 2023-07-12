#ifndef __eserc1classes__
#define __eserc1classes__

#include "random.h"
#include "eserc1functions.h"

using namespace std;

class r {

    public:

        double Eval(double x){
            return x;
        }
        double Evalerror(double x){
            return pow(x - 0.5, 2.);
        }

};

class BlockAverage{

    public:

        BlockAverage(r& fun, int n = 2, int b = 20):
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

        void Measure(double x){

            walker[iav] = function.Eval(x);
            walker[ierr] = function.Evalerror(x);
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block
        
            ofstream Average, Error; 
            const int wd=12;
            
            cout << "Block number " << iblk << endl;

            Average.open("output_average.dat",ios::app);
            Error.open("output_error.dat", ios::app);

            stima = blk_av[iav]/blk_norm; // f(r) = r
            glob_av[iav] += stima;
            glob_av2[iav] += stima*stima;
            errore = error(glob_av[iav],glob_av2[iav],iblk);

            Average << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iav]/(double)iblk << setw(wd) << errore << endl;

            stima = blk_av[ierr]/blk_norm; // f(r) = (r - 0.5)^2
            glob_av[ierr] += stima;
            glob_av2[ierr] += stima*stima;
            errore = error(glob_av[ierr],glob_av2[ierr],iblk);

            Error << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[ierr]/(double)iblk << setw(wd) << errore << endl;

            cout << "----------------------------" << endl << endl;

            Average.close();
            Error.close();

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

        r function;

};


#endif 