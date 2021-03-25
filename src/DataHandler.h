 #ifndef  __DATAHANDLER_H_
#define   __DATAHANDLER_H_

#include <string.h>
#include<iostream>

// #include "Timer.h"
#include "Supporter.h"

class pMInstance{
public:
	char 													_czfile_name[200];
	char 													_facilityfile_name[200];
	char														_transcostfile_name[200];
	int 														_nb_cz_nodes;
	int 														_nb_facility_nodes;

	double *												_cz_demands = nullptr;
	/* File type I: given coordinates (x,y) in 2D space */
	int														_nb_dims;
	double ** 											_cz_coors = nullptr;
	double ** 											_facility_coors = nullptr;
	/* File type II: given cost table and demand list (no coordinates)*/
	double **												_transcost_matx = nullptr;
	/* File type III: given latitude and longtitude in geographical space */
	double **												_cz_latlong_coors = nullptr;
	double ** 											_facility_latlong_coors = nullptr;

	/*read file <cz or facility>_<format>
									[X: coordinate x]
									[Y: corrdinate y]
			[cz]					[D: demand value]
			[facility]			[Table: cost table]
									[La: latitude]
									[Lo: longtitude]                     */
	void read_czfile_XYD(const char *filename_);
	void read_facilityfile_XY(const char *filename_);
	void read_transcost_Table(const char *filename_, const char *filename2_);
	void read_demands_D(const char * filename_);
	void read_czfile_LaLoD(const char*filename_);
	void read_facilityfile_LaLo(const char*filename_);

	virtual ~pMInstance(void);
	inline void error(const char* s1, const char *s2 = ""){
		printf("\nERROR:\n%s %s\n", s1, s2);
		exit(1);
	};

	inline long long atoll_(const char* str){
	    long long i = 0;
	    while(strcmp(str, "") != 0){
	       i*=10;
	       i += *str++ - '0';
	    }
	    return i;
	}
};


#endif  // __DATAHANDLER_H_