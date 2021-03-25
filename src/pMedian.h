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




#define EPSILON   (float)0.000001
#define INF std::numeric_limits<float>::infinity() ;

/*Here,  */
typedef tuple<float, vector<float>, vector<int>> TYPE_ZRY;
// typedef pair<float, vector<vector<int>>> TYPE_ZX;



struct Arc{
		float 											_dist;
		int 												_idx;
};


class pMedian{
	public:
		pMInstance	*								_inst = nullptr;
		int												_p;
		int 												_nb_customers;
		int 												_nb_facilities;

		float **										_dist_matx = nullptr;
		float **										_w_dist_matx = nullptr; // weight distance matrix: demand[j]* dist[i,j]
		Arc**											_wdist_matx_sorted = nullptr;	


		void initialize_by_fmt_XY(pMInstance & inst_, int p_);
		void initialize_by_fmt_Table(pMInstance & inst_, int p_);
		void initialize_by_fmt_LaLo(pMInstance & inst_, int p_);


		virtual ~pMedian(void);

		TYPE_ZRY calc_lower_bound(const vector<float> & lambda_);
		double calc_upper_bound(const vector<int> &Y_);
		void Lagrangian_algo();
		void initialize_LagMultipliers(vector<float> & lambda_);
		void stabilize_LagMultipliers(vector<float> & lambda_, int itr);

		void calc_weighted_dist_matx(float (*calc_dist)(pMInstance*, int, int));
		static float calc_eucl_dist(pMInstance*, int fac_idx_, int cz_idx_);
		static float to_Radians(const float degree);
		static float calc_geo_dist(pMInstance*, int fac_idx_, int cz_idx_);

};

#endif  //  __PMEDIAN_H_
