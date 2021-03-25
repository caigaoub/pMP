#include <string>
#include <fstream>

#include "DataHandler.h"
using namespace std;

/* Read customer location & demand file */
void pMInstance::read_czfile_XYD(const char *filename_){
	strcpy(_czfile_name, filename_);
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open customer file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading customer file ...");
	}
	infile >> _nb_cz_nodes >> _nb_dims;
	// Format::pretty_print( "Nodes: " + to_string(_nb_cz_nodes) + ", Dimensions: " + to_string(_nb_dims)) ;
	// read node coordinates and demands 
	_cz_coors = new double*[_nb_cz_nodes];
	for(int i = 0; i < _nb_cz_nodes; i++){
		_cz_coors[i] = new double[_nb_dims];
	}
	_cz_demands = new double[_nb_cz_nodes];
	for(int i = 0; i < _nb_cz_nodes; i++) {
		infile >> _cz_coors[i][0] >> _cz_coors[i][1] >> _cz_demands[i];
	}
	infile.close();

	Format::pretty_print("The first 10 rows are: " );
	if(1){
		for(int i=0; i< min(_nb_cz_nodes,10); i++) {
			printf("%5d", i+1);
			for(int j=0; j< _nb_dims; j++)
				printf(" %10.2f", _cz_coors[i][j]);
			printf(" %10.2f", _cz_demands[i]);
			printf("\n");
		}
	}
	Format::pretty_print("Customer-info file reading is done in " + to_string(Timer::CPU_time() - tt) + " seconds");

}


/* Read facilty location file */
void pMInstance::read_facilityfile_XY(const char *filename_){
	strcpy(_facilityfile_name, filename_);
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open facility file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading facility file ...");
	}
	infile >> _nb_facility_nodes >> _nb_dims;
	// Format::pretty_print( "Nodes: " + to_string(_nb_cz_nodes) + ", Dimensions: " + to_string(_nb_dims)) ;
	// read node coordinates and demands 
	_facility_coors = new double*[_nb_facility_nodes];
	for(int i = 0; i < _nb_facility_nodes; i++){
		_facility_coors[i] = new double[_nb_dims];
	}
	for(int i = 0; i < _nb_facility_nodes; i++) {
		infile >> _facility_coors[i][0] >> _facility_coors[i][1];
	}
	infile.close();

	Format::pretty_print("The first 10 rows are: " );
	if(1){
		for(int i=0; i< min(_nb_facility_nodes,10); i++) {
			printf("%5d", i+1);
			for(int j=0; j< _nb_dims; j++)
				printf(" %10.2f", _facility_coors[i][j]);
			printf("\n");
		}
	}
	Format::pretty_print("Facility-info file reading is done in " + to_string(Timer::CPU_time() - tt) + " seconds");

}


void pMInstance::read_transcost_Table(const char * filename_, const char * filename2_){
	strcpy(_transcostfile_name, filename_);
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open trans-cost file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading trans-cost file ...");
	}
	infile >> _nb_facility_nodes >> _nb_cz_nodes;
	// read node coordinates and demands 
	_transcost_matx = new double*[_nb_facility_nodes];
	for(int i = 0; i < _nb_facility_nodes; i++){
		_transcost_matx[i] = new double[_nb_cz_nodes];
	}
	for(int i = 0; i < _nb_facility_nodes; i++) {
		for(int j=0; j<_nb_cz_nodes; j++){
			infile >> _transcost_matx[i][j];
		}
	}
	infile.close();

	Format::pretty_print("The first 10 rows are: " );
	if(1){
		for(int i=0; i< min(_nb_facility_nodes,10); i++) {
			printf("%5d", i+1);
			for(int j=0; j< _nb_cz_nodes; j++)
				printf(" %5.2f", _transcost_matx[i][j]);
			printf("\n");
		}
	}
	Format::pretty_print("Facility-info file reading is done in " + to_string(Timer::CPU_time() - tt) + " seconds");
	read_demands_D(filename2_);
}

void pMInstance::read_demands_D(const char * filename_){
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open demand file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading demand file ...");
	}
	// readdemands 
	_cz_demands = new double[_nb_cz_nodes];
	for(int j=0; j<_nb_cz_nodes; j++){
		infile >> _cz_demands[j];
	}
	infile.close();
}

void pMInstance::read_czfile_LaLoD(const char*filename_){
	strcpy(_czfile_name, filename_);
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open customer file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading customer file ...");
	}
	infile >> _nb_cz_nodes >> _nb_dims;
	// Format::pretty_print( "Nodes: " + to_string(_nb_cz_nodes) + ", Dimensions: " + to_string(_nb_dims)) ;
	// read node coordinates and demands 
	_cz_latlong_coors = new double*[_nb_cz_nodes];
	for(int i = 0; i < _nb_cz_nodes; i++){
		_cz_latlong_coors[i] = new double[_nb_dims];
	}
	_cz_demands = new double[_nb_cz_nodes];
	for(int i = 0; i < _nb_cz_nodes; i++) {
		infile >> _cz_latlong_coors[i][0] >> _cz_latlong_coors[i][1] >> _cz_demands[i];
	}
	// for(int i = 0; i < _nb_cz_nodes; i++) {
	// 	_cz_demands[i]/=100000.0;
	// }
	infile.close();

	Format::pretty_print("The first 10 rows are: " );
	if(1){
		for(int i=0; i< min(_nb_cz_nodes,10); i++) {
			printf("%5d", i+1);
			for(int j=0; j< _nb_dims; j++)
				printf(" %10.2f", _cz_latlong_coors[i][j]);
			printf(" %10.2f", _cz_demands[i]);
			printf("\n");
		}
	}
	Format::pretty_print("Customer-info file reading is done in " + to_string(Timer::CPU_time() - tt) + " seconds");

}



/* Read facilty location file */
void pMInstance::read_facilityfile_LaLo(const char *filename_){
	strcpy(_facilityfile_name, filename_);
	double tt = Timer::CPU_time();
	// open instance file
	ifstream infile(filename_, ifstream::in);
	if(!infile.is_open()){
		Format::pretty_print("Fail to open facility file\n");
		exit(0);
	}else{
		Format::pretty_print("Start reading facility file ...");
	}
	infile >> _nb_facility_nodes >> _nb_dims;
	// Format::pretty_print( "Nodes: " + to_string(_nb_cz_nodes) + ", Dimensions: " + to_string(_nb_dims)) ;
	// read node coordinates and demands 
	_facility_latlong_coors = new double*[_nb_facility_nodes];
	for(int i = 0; i < _nb_facility_nodes; i++){
		_facility_latlong_coors[i] = new double[_nb_dims];
	}
	for(int i = 0; i < _nb_facility_nodes; i++) {
		infile >> _facility_latlong_coors[i][0] >> _facility_latlong_coors[i][1];
	}
	infile.close();

	Format::pretty_print("The first 10 rows are: " );
	if(1){
		for(int i=0; i< min(_nb_facility_nodes,10); i++) {
			printf("%5d", i+1);
			for(int j=0; j< _nb_dims; j++)
				printf(" %10.2f", _facility_latlong_coors[i][j]);
			printf("\n");
		}
	}
	Format::pretty_print("Facility-info file reading is done in " + to_string(Timer::CPU_time() - tt) + " seconds");
}



pMInstance::~pMInstance(void){
	if(_cz_coors){ // free customer coords & demands
		for(int i=0; i<_nb_cz_nodes; i++){
			delete[] _cz_coors[i];
		}
		delete[] _cz_coors;
		delete[] _cz_demands;
	}
	if(_facility_coors){ // free facility coords
		for(int i=0; i<_nb_facility_nodes; i++){
			delete[] _facility_coors[i];
		}
		delete[] _facility_coors;
	}

	if(_transcost_matx){ // free trans-cost table & customer demand array
		for(int i=0; i<_nb_facility_nodes; i++){
			delete[] _transcost_matx[i];
		}
		delete[] _transcost_matx;
		delete[] _cz_demands;
	}
	if(_cz_latlong_coors){ // free customer coords & demands
		for(int i=0; i<_nb_cz_nodes; i++){
			delete[] _cz_latlong_coors[i];
		}
		delete[] _cz_latlong_coors;
		delete[] _cz_demands;
	}
	if(_facility_latlong_coors){ // free facility coords
		for(int i=0; i<_nb_facility_nodes; i++){
			delete[] _facility_latlong_coors[i];
		}
		delete[] _facility_latlong_coors;
	}
}