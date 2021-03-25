#include "pMedian.h"


/**********************************************************
    ****************  MIP formluation    ************************

				MIN \sum_I \sum_J {d_j*c_{ij}*x_{ij}}    (1)
			s.t. 
				\sum_{i\in I} x_{ij} = 1, \forall j \in J,      (2)
				\sum_{j\in J} y_i = p,                                  (3)
				x_{ij} <= y_i, \forall i \in I, j \in J,             (4)
				y_i binary,  \forall i\in I,                               (5)
				x_{ij} >= 0, \forall i\in I, j\in J                   (6)

***********************************************************
***********************************************************/


void pMedian::initialize_by_fmt_XY(pMInstance &inst_, int p_){
	// read instance parameters
	this->_inst = &inst_;
	this->_nb_customers = _inst->_nb_cz_nodes;
	this->_nb_facilities = _inst->_nb_facility_nodes;
	this->_p = p_;
	// calculate distance matrx (include weighted distance matrx)
	calc_weighted_dist_matx(pMedian::calc_eucl_dist);
}

void pMedian::initialize_by_fmt_Table(pMInstance &inst_, int p_){
	// read instance parameters
	this->_inst = &inst_;
	this->_nb_customers = _inst->_nb_cz_nodes;
	this->_nb_facilities = _inst->_nb_facility_nodes;
	this->_p = p_;
	// calculate distance matrx (include weighted distance matrx)
	// calc_weighted_dist_matx();
	this->_dist_matx = _inst->_transcost_matx;
	_w_dist_matx = new float*[_nb_facilities];
	for(int i=0; i<_nb_facilities; i++){
		_w_dist_matx[i] = new float[_nb_customers];
	}
	/* sorted wdist matrix*/
	_wdist_matx_sorted = new Arc*[_nb_customers];
	for(int j=0; j<_nb_customers; j++){
		_wdist_matx_sorted[j] = new Arc[_nb_facilities];
	}
	for(int i=0; i<_nb_facilities; i++){
		for(int j=0; j <_nb_customers; j++){
			_w_dist_matx[i][j] =  _inst->_cz_demands[j] * _inst->_transcost_matx[i][j];
			_wdist_matx_sorted[j][i] = {_w_dist_matx[i][j], i};
		}
	}
	// sort the matrix
	for(int j=0; j<_nb_customers; j++){
		std::sort(_wdist_matx_sorted[j], _wdist_matx_sorted[j]+_nb_facilities,[](const Arc &a,  const Arc &b){return a._dist < b._dist;});
	}

}

void pMedian::initialize_by_fmt_LaLo(pMInstance & inst_, int p_){
	this->_inst = &inst_;
	this->_nb_customers = _inst->_nb_cz_nodes;
	this->_nb_facilities = _inst->_nb_facility_nodes;
	this->_p = p_;
	calc_weighted_dist_matx(pMedian::calc_geo_dist);
}


/* pre-calculate the weighted distance matrix */
void pMedian::calc_weighted_dist_matx(float (*calc_dist)(pMInstance*, int, int)){
	_w_dist_matx = new float*[_nb_facilities];
	// _dist_matx = new float*[_nb_facilities];
	for(int i=0; i<_nb_facilities; i++){
		_w_dist_matx[i] = new float[_nb_customers];
		// _dist_matx[i] = new float[_nb_customers];
	}
	/* sorted wdist matrix*/
	_wdist_matx_sorted = new Arc*[_nb_customers];
	for(int j=0; j<_nb_customers; j++){
		_wdist_matx_sorted[j] = new Arc[_nb_facilities];
	}
	for(int i=0; i<_nb_facilities; i++){
		for(int j=0; j <_nb_customers; j++){
			// _dist_matx[i][j] = calc_eucl_dist(i, j); 
			_w_dist_matx[i][j] = _inst->_cz_demands[j] * calc_dist(_inst, i, j);
			_wdist_matx_sorted[j][i] = {_w_dist_matx[i][j], i};
		}
	}
		// sort the matrix
	for(int j=0; j<_nb_customers; j++){
		std::sort(_wdist_matx_sorted[j], _wdist_matx_sorted[j]+_nb_facilities,[](const Arc &a,  const Arc &b){return a._dist < b._dist;});
	}


	if(0){
		for(int i=0; i<_nb_facilities; i++){
			for(int j=0; j <_nb_customers; j++){
					printf(" %5.2f", _w_dist_matx[i][j]); 
			}
			printf("\n"); 
		}
	}
}

/*  euclidean distance function */
float pMedian::calc_eucl_dist(pMInstance * inst_, int fac_idx_, int cz_idx_){
	float dist = pow(inst_->_facility_coors[fac_idx_][0] - inst_->_cz_coors[cz_idx_][0], 2);
	dist += pow(inst_->_facility_coors[fac_idx_][1] - inst_->_cz_coors[cz_idx_][1], 2);
	dist = max(dist, EPSILON); // avoid numerical error when custmer and facility are in the same location 
	return pow(dist,0.5);
}


/*  geographical distance function */
float pMedian::calc_geo_dist(pMInstance * inst_,  int fac_idx_, int cz_idx_){
    float  lat1 = pMedian::to_Radians(inst_->_facility_latlong_coors[fac_idx_][0]);
    float  long1 = pMedian::to_Radians(inst_->_facility_latlong_coors[fac_idx_][1]);
    float lat2 = pMedian::to_Radians(inst_->_cz_latlong_coors[cz_idx_][0]);
    float long2 = pMedian::to_Radians(inst_->_cz_latlong_coors[cz_idx_][1]);
      
    // Haversine Formula
    float dlong = long2 - long1;
    float dlat = lat2 - lat1;
  
    float ans = pow(sin(dlat / 2), 2) + cos(lat1)*cos(lat2)* pow(sin(dlong / 2), 2);
  
    ans = 2 * asin(sqrt(ans));
  
    // Radius of Earth in 
    // Kilometers, R = 6371
    // Use R = 3956 for miles
    float R = 3956.0;
    // Calculate the result
    ans = ans * R;
    return ans;
}

float pMedian::to_Radians(const float degree){
    float one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/

void pMedian::Lagrangian_algo(){
	vector<float> lambda(_nb_customers); // initialize lagr. multipliers
	initialize_LagMultipliers(lambda);
	vector<vector<int>> best_X; // assignment variables of MIP: X
	vector<int> best_Y; // indicator variables that decide which subset of facility to open: Y
	float Z_UB = INF; // upper bound
	float Z_LB = -INF; // lower bound
	float theta = 1.;
	int MAX_ITER = 100000;
	int itr = 1;
	int lb_not_updated = 0;
	while(true){
		// obtaining the lower and upper bounds
		auto ret = calc_lower_bound(lambda);
		float Z_D = get<0>(ret);
		vector<float> R = get<1>(ret);

		auto ret2 = calc_upper_bound(get<2>(ret));
		float Z = ret2.first;
		// update the upper bound
		if(Z < Z_UB){
			Z_UB = Z;
			best_X = ret2.second;
			best_Y = get<2>(ret);
		}
		//update the lower bound
		if(Z_D > Z_LB){
			Z_LB = Z_D;
			lb_not_updated = 0;
		}else{
			lb_not_updated++;
			if (lb_not_updated >= 30){
				theta /= 2.0;
				lb_not_updated = 0;
			}
		}
		if(theta < 0.0000005){
			
			Format::pretty_print("Iteration stops due to theta = 0 \n");
			break;
		}
		// determine the step size and update the multiplier
		// vector<float> residual(_nb_customers, 0.0);
		// for(int j=0; j<_nb_customers; j++){
		// 	residual[j] = 1.0;
		// 	for(int i=0; i <_nb_facilities; i++){
		// 		residual[j] -= X_D[i][j];
		// 	}
		// }
		float residual_sq_sum = 0.0;
		for(int j=0; j <_nb_customers; j++){
			residual_sq_sum += pow(R[j],2);
		}
		float t = theta * (Z_UB - Z_D)/ residual_sq_sum;
		// cout << "theta = " << theta <<endl;
		// update lambda
		for(int j=0; j<_nb_customers; j++){
			lambda[j] += t * R[j];	
		}
		// stabilize_LagMultipliers(lambda, itr+2);
		itr++;
		// compute the optimality gap
		float opt_gap = (Z_UB - Z_LB)/Z_LB;
		if (opt_gap < EPSILON){
			Format::pretty_print("optimality gap is met " + to_string(opt_gap *100) + "%");
			Format::pretty_print(to_string(Z_LB) + ", "+ to_string(Z_UB) + ", " + to_string(opt_gap*100) + "%");

			break;
		}
		// cout << Z_LB << ", " << Z_UB << ", " << opt_gap*100 << "%" << endl;
		Format::pretty_print(to_string(Z_LB) + ", "+ to_string(Z_UB) + ", " + to_string(opt_gap*100) + "%");

	}


}

void pMedian::initialize_LagMultipliers(vector<float> & lambda_){
	for(int j=0; j<_nb_customers; j++){
		lambda_[j] = _wdist_matx_sorted[j][0]._dist;
	}
}

void pMedian::stabilize_LagMultipliers(vector<float> & lambda_, int itr){
	// vector<pair<float,int>> lambda_ub(_nb_customers);
	int idx = min(itr, _nb_facilities-1);
	for(int j=0; j <_nb_customers; j++){
		lambda_[j] = min(lambda_[j], _wdist_matx_sorted[idx][j]._dist);
	}
}

/*
Given a set of multipliers, optimality of solving the lag relaxation 
will give the "best" lower bound on MIP
	*/
TYPE_ZRY pMedian::calc_lower_bound(const vector<float> & lambda_){
	// obtain v
	vector<pair<float, int>> v(_nb_facilities);
	for(int i=0; i<_nb_facilities; i++){
		v[i] = make_pair(0.0, i);
	}
	for(int i=0; i < _nb_facilities; i++){
		for(int j=0; j< _nb_customers; j++){
			v[i].first += min((float)0.0, _w_dist_matx[i][j] - lambda_[j]); 
		}
	}

	std::sort(v.begin(), v.end(),[](const pair<float,int> &a,  const pair<float,int> &b){return a.first < b.first;});

	// determine Y
	vector<int> Y(_nb_facilities, 0);
	for(int i=0; i< this->_p; i++){
		Y[v[i].second] = 1;
	}
	// lastest way
	// vector<vector<int>> X(_nb_facilities , vector<int> (_nb_customers, 0));
	Arc tmpE;
	float Z_D = 0.0;
	vector<float> R(_nb_customers, 1.0);
	for(int j=0; j< _nb_customers; j++){
		Z_D += lambda_[j];
		for(int k=0; k < _nb_facilities; k++){
			tmpE = _wdist_matx_sorted[j][k];
			if( Y[tmpE._idx] ==1){
				if(tmpE._dist - lambda_[j] < 0){
					// X[tmpE._idx][j] = 1;
					R[j] -= 1.0;
					Z_D += tmpE._dist - lambda_[j];	
				}else{
					break;
				}			
			}
		}
	}

	// New way: determine X
	// vector<vector<int>> X(_nb_facilities , vector<int> (_nb_customers, 0));
	// for(int k=0; k < this->_p; k++){
	// 	for(int j=0; j< _nb_customers; j++){
	// 		if( _w_dist_matx[v[k].second][j] - lambda_[j] < 0){
	// 			X[v[k].second][j] = 1;
	// 		}
	// 	}
	// }
	// float Z_D = 0.0;
	// for(int j=0; j <_nb_customers; j++){
	// 	Z_D += lambda_[j];
	// }
	// for(int k=0; k <this->_p; k++){
	// 	for(int j=0; j<_nb_customers; j++){
	// 		if(X[v[k].second][j] == 1)
	// 			Z_D += _w_dist_matx[v[k].second][j] -lambda_[j];
	// 	}
	// }

	// Old way: determine X
	// vector<vector<int>> X(_nb_facilities , vector<int> (_nb_customers, 0));
	// for(int i=0; i < _nb_facilities; i++){
	// 	for(int j=0; j< _nb_customers; j++){
	// 		if( Y[i]==1 && _w_dist_matx[i][j] - lambda_[j] < 0){
	// 			X[i][j] = 1;
	// 		}
	// 	}
	// }
	// // compute Z_D (the objective function of the lagr relaxation)
	// float Z_D = 0.0;
	// for(int j=0; j <_nb_customers; j++){
	// 	Z_D += lambda_[j];
	// 	for(int i=0; i<_nb_facilities; i++){
	// 		Z_D += _w_dist_matx[i][j]*X[i][j] -lambda_[j]*X[i][j];
	// 	}
	// }
	return make_tuple(Z_D, R, Y);
}

TYPE_ZX pMedian::calc_upper_bound(const vector<int> &Y_){
	//given Y, compute X
	vector<vector<int>> X(_nb_facilities , vector<int> (_nb_customers, 0));
	float tmp = INF;
	int nearest_fac_idx = -1;
	for(int j=0; j<_nb_customers; j++){
		 tmp = INF;
		 nearest_fac_idx = -1;
		for(int i=0; i <_nb_facilities; i++){
			if(Y_[i] == 1 && _w_dist_matx[i][j] < tmp){
				tmp = _w_dist_matx[i][j];
				nearest_fac_idx = i;
			}
		}
		X[nearest_fac_idx][j] = 1;
	}
	// compute Z
	float Z = 0.0;
	for(int i=0; i<_nb_facilities; i++){
		for(int j=0; j<_nb_customers; j++){
			Z +=_w_dist_matx[i][j] * X[i][j];
		}
	}
	// Format::pretty_print("current feasible solution: " + to_string(Z));
	return	make_pair(Z, X);
}



pMedian::~pMedian(void){
	if(_w_dist_matx){ // free trans-cost table & customer demand array
		for(int i=0; i<_nb_facilities; i++){
			delete[] _w_dist_matx[i];
		}
		delete[] _w_dist_matx;
		// delete[] _cz_demands;
	}
	if(_wdist_matx_sorted){ 
		for(int i=0; i<_nb_customers; i++){
			delete[] _wdist_matx_sorted[i];
		}
		delete[] _wdist_matx_sorted;
	}
	// if(!_dist_matx){ // 
	// 	for(int i=0; i<_nb_facilities; i++){
	// 		delete[] _dist_matx[i];
	// 	}
	// 	delete[] _dist_matx;
	// }
}