//First example

// #include "mpi.h"//"/opt/homebrew/Cellar/mpich/4.1.1_1/include/mpi.h" //include "mpi.h"

// #include <iostream>
// using namespace std;
// int main(int argc, char* argv[])
// {
//     int size, rank; 
//     MPI_Init(&argc,&argv); 
    
//     MPI_Comm_size(MPI_COMM_WORLD, &size); 
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//     cout<<" Sono il nodo "<<rank<<" dei "<<size<<" che hai utilizzato!"<<endl;
//     MPI_Finalize();
//     return 0; 
// }

//############################

//MPI_BCAST

// #include "mpi.h" 
// #include <iostream> 

// using namespace std;
// int main(int argc, char* argv[])
// {
// int size, rank; 
// MPI_Init(&argc,&argv); 
// MPI_Comm_size(MPI_COMM_WORLD, &size); 
// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// int my_values[3]; 
// for(int i=0;i<3;i++){
// if(rank==0) my_values[i]=i+1;
// else my_values[i]=0;
// }
// cout<< "Prima: "<< my_values[0]<< " "<< my_values[1]<< " "<< my_values[2]<< " per il processo "<< rank<< endl;
// MPI_Bcast(my_values,3,MPI_INTEGER,0, MPI_COMM_WORLD); 
// cout<< "Dopo: "<< my_values[0]<< " "<< my_values[1]<< " "<< my_values[2]<< " per il processo "<< rank<< endl;
// MPI_Finalize();
// return 0; 
// }

//###################################

//MPI_GATHER

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// int main(int argc, char* argv[]) {
// int size, rank;
// MPI_Init(&argc,&argv); 
// MPI_Comm_size(MPI_COMM_WORLD, &size); 
// MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
// if(size>3){cout<<"Hai scelto troppi processi"<<endl;
// return 1;} 
// int irecv[3];
// for(int i=0;i<3;i++) {
// irecv[i]=0; }
// int isend = rank + 1;

// MPI_Gather(&isend,1,MPI_INTEGER,irecv,1,MPI_INTEGER,0, MPI_COMM_WORLD);
// if(rank==0) cout<< "irecv: " <<irecv[0] <<" " <<irecv[1] <<" " <<irecv[2] <<endl;
// MPI_Finalize();
// return 0; }


//########################################


//MPI_REDUCE

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// int main(int argc, char* argv[]) {
// int size, rank; 
// MPI_Init(&argc,&argv); 
// MPI_Comm_size(MPI_COMM_WORLD, &size); 
// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// int isend[2],irecv[2];
// for(int i=0;i<2;i++) isend[i]=rank+i+1;

// MPI_Reduce(&isend[0],&irecv[0],1,MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD);
// MPI_Reduce(&isend[1],&irecv[1],1,MPI_INTEGER, MPI_PROD,0,MPI_COMM_WORLD);

// if(rank==0)
// cout<<"irecv: "<<irecv[0]<<" "<<irecv[1]<<endl;

// MPI_Finalize();
// return 0; }

//#######################################


//MPI_COMM_SPLIT

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// int main(int argc, char* argv[]) {
// int size, rank;
// MPI_Init(&argc,&argv);
// MPI_Comm_size(MPI_COMM_WORLD, &size); 
// MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
// if(size!=4){cout<<"Servono 4 processi, non "<<size<<"!!"<<endl; return 1;}
// int icolor, ikey;
// if(rank==0){icolor=1;ikey=2;} 
// if(rank==1){icolor=1;ikey=1;} 
// if(rank==2){icolor=2;ikey=1;} 
// if(rank==3){icolor=2;ikey=2;}
// MPI_Comm nuovocom; 
// MPI_Comm_split(MPI_COMM_WORLD,icolor,ikey,&nuovocom); 
// int newsize,newrank;
// MPI_Comm_size(nuovocom, &newsize); 
// MPI_Comm_rank(nuovocom, &newrank);
// cout<<"Ero: "<<rank<<" di "<<size<<" ... e adesso sono: "<<newrank<<" di "<<newsize<<endl;
// MPI_Finalize();
// return 0;
// }


//#####################################

//MPI_SEND/MPI_RECV

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// int main(int argc, char* argv[]) {
// int size, rank; MPI_Init(&argc,&argv); MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Status stat;
// int itag=1;
// int imesg = rank; if(rank==1)
// MPI_Send(&imesg,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD); else if(rank==0)
// MPI_Recv(&imesg,1,MPI_INTEGER,1,itag,MPI_COMM_WORLD, &stat);
// cout<<"messaggio = "<<imesg<<endl;
// MPI_Finalize(); return 0;
// }



//Bidirectional communications

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// const int n = 100; // try to increase n int main(int argc, char* argv[]){
// int size, rank;
// MPI_Init(&argc,&argv);
// MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// MPI_Status stat1, stat2;
// int* imesg = new int[n]; int* imesg2 = new int[n]; int itag=1; int itag2=2;
// for(int i=0;i<n;i++){imesg[i]=rank; imesg2[i]=rank+1;} if(rank==1){MPI_Send(&imesg[0],n,
// MPI_INTEGER,0,itag,MPI_COMM_WORLD); MPI_Recv(&imesg2[0],n,
// MPI_INTEGER,0,itag2, MPI_COMM_WORLD,&stat2); cout<<"messaggio = "<<imesg2[0]<<endl;}
// else if(rank==0){MPI_Send(&imesg2[0],n, MPI_INTEGER,1,itag2, MPI_COMM_WORLD);
// MPI_Finalize(); return 0;}



//Bidirectional communications #2

// #include "mpi.h"
// #include <iostream>
// using namespace std;
// const int n = 100; // try to increase n int main(int argc, char* argv[]){
// int size, rank; MPI_Init(&argc,&argv); MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Status stat1, stat2;
// MPI_Request req;
// int* imesg = new int[n]; int* imesg2 = new int[n]; int itag=1; int itag2=2;
// for(int i=0;i<n;i++){imesg[i]=rank; imesg2[i]=rank+1;} if(rank==1){MPI_Isend(&imesg[0],n,
// MPI_INTEGER,0,itag, MPI_COMM_WORLD,&req); MPI_Recv(&imesg2[0],n,MPI_INTEGER,0,itag2, MPI_COMM_WORLD,&stat2);
// cout<<"messaggio = "<<imesg2[0]<<endl;}
// else if(rank==0){MPI_Send(&imesg2[0],n, MPI_INTEGER,1,itag2,MPI_COMM_WORLD);
// MPI_Finalize(); return 0;}



//MPI_WTIME

// #include "mpi.h" #include <iostream>
// using namespace std;
// int main(int argc, char* argv[]){ int size, rank; MPI_Init(&argc,&argv); MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Status stat;
// double tstart = MPI_Wtime(); int n = 100;
// int* imesg = new int[n]; int sum=0;
// for(int i=0;i<n;i++){
// imesg[i]=rank;
// if(rank==1) MPI_Send(&imesg[0],n,
// MPI::INTEGER,0,i,MPI_COMM_WORLD); else if(rank==0) MPI_Recv(&imesg[0],n,MPI_INTEGER,
// 1,i,MPI_COMM_WORLD,&stat);
// sum += imesg[i];} double tend = MPI_Wtime();
// double dt = tend - tstart;
// cout<<"io sono "<<rank<<"; somma = "<<sum<<"; tempo = "<<dt<<endl;
// MPI_Finalize();
// return 0;}