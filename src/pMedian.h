 #ifndef  __PMEDIAN_H_
#define  __PMEDIAN_H_


#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <limits>
#include <functional>
#include <bits/stdc++.h>

#include "DataHandler.h"




#define EPSILON   0.000001
#define INF std::numeric_limits<double>::infinity() ;

/*Here, Z, X, Y are variables we define for the MIP */
typedef tuple<double, vector<vector<int>>, vector<int>> TYPE_ZXY;
typedef pair<double, vector<vector<int>>> TYPE_ZX;



struct Arc{
		double 											_dist;
		int 												_idx;
};


class pMedian{
	public:
		pMInstance	*								_inst = nullptr;
		int												_p;
		int 												_nb_customers;
		int 												_nb_facilities;

		double **										_dist_matx = nullptr;
		double **										_w_dist_matx = nullptr; // weight distance matrix: demand[j]* dist[i,j]
		Arc**											_wdist_matx_sorted = nullptr;	


		void initialize_by_fmt_XY(pMInstance & inst_, int p_);
		void initialize_by_fmt_Table(pMInstance & inst_, int p_);
		void initialize_by_fmt_LaLo(pMInstance & inst_, int p_);


		virtual ~pMedian(void);

		TYPE_ZXY calc_lower_bound(const vector<double> & lambda_);
		TYPE_ZX calc_upper_bound(const vector<int> &Y_);
		void Lagrangian_algo();
		// static bool sort_by_sec(const pair<double,int> &a,  const pair<double,int> &b);

		void calc_weighted_dist_matx(double (*calc_dist)(pMInstance*, int, int));
		static double calc_eucl_dist(pMInstance*, int fac_idx_, int cz_idx_);
		static double to_Radians(const double degree);
		static double calc_geo_dist(pMInstance*, int fac_idx_, int cz_idx_);

};

#endif  //  __PMEDIAN_H_
