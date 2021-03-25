
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include<iostream>


#include "DataHandler.h"
#include "pMedian.h"
#include "Supporter.h"

using namespace std;

int main(int argc, char* argv[]){
switch(argv[1][0])
	{
	case 'V': // cheching solution
		{
			break;
		}
	case 'L': // only subgradient
		{ 
			pMInstance instance;
			
			double tt = Timer::CPU_time();

			try{
				instance.read_czfile_LaLoD(argv[2]); // read instance
				instance.read_facilityfile_LaLo(argv[3]); // read instance
				pMedian	mdl;
				int p_ = atoi(argv[4]);
				mdl.initialize_by_fmt_LaLo(instance, p_);
				
				// instance.read_transcost_Table(argv[2], argv[3]);
				// pMedian	mdl;
				// int p_ = atoi(argv[4]);
			 // 	mdl.initialize_by_fmt_Table(instance, p_);



				mdl.Lagrangian_algo();
				// pmed.p = atoi(argv[3]); // parameter p 
				// pmed.bub = atof(argv[4]); // best upper bound
				// pmed.r = atoi(argv[5]); // parametrer r
				// pmed.MEMSIZE = atoll_(argv[6]); // memery size setting 
				Format::pretty_print("Instance is solved in " + to_string(Timer::CPU_time() - tt) + " seconds");

			}
			catch(...){
				cout << "Wrong command arguments." << endl;
				exit(0);
			}
			// pmed.lang_set1();
			// pmed.PRIMALFREQ = 0;
			// pmed.writedist();
			// //pmed.readdist();
			// double tt1 = CPU_time();
			// pmed.lang_init();
			// pmed.lang();
			// pmed.lang_final();

			// printf("Total time %.2f\n", CPU_time() - tt);
			// printf("Lag time %.2f\n", CPU_time() - tt1);
			// FILE *ff = fopen("res.txt", "a");
			// fprintf(ff, "%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.4f\t%.2f%.2f\n", pmed.probname, 
			// 	pmed.m, pmed.n, pmed.p, pmed.blb, pmed.bub, (pmed.bub-pmed.blb)/pmed.bub * 100., CPU_time() - tt,CPU_time() - tt1);
			// fclose(ff);
			break;
		}
	default:
		{
			cout << "Unknown method lable." << endl;
			exit(0);
		}
	}
	return 0;

}