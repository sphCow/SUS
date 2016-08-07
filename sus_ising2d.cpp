#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <string>

using namespace std;

enum Direction {L=0,R=1};

//***** GLOBAL PARAMETERS *****//
int nx;
int N;
double T;
long mcs;
vector<int> lattice; //this contains spin values, +/- 1 int ix = index%nx;int iy = index/nx;

//rng
std::random_device rd;
std::mt19937 rng(123456);
std::uniform_real_distribution<double> dist(0.0,1.0);
//rng; //std::random_device()());

//files
ofstream flat;
ofstream fthermo;
ofstream fthermo2;

//***** Generate lattice ******//
void generate_lattice(int init) {
    /* init =  1 => all lattice +1
     *      = -1 => all lattice -1
     *      =  0 => random  */
    
    for(int i=0; i<N; i++) {
        if(init == 1 || init == -1) lattice[i] = init;
        else {
            lattice[i] = ( dist(rng)>0.5 ? 1 : -1);
        }
        
    }
}

//******* print lattice ******//
void print_lattice() {
    for(int i=0; i<N; i++) {
        flat << i%nx << " " << int(i/nx) << " " 
             << (lattice[i]==1 ? "+" : "-") << endl;
    }
    flat << endl << endl;
}

//******* calculate energy ******//
int get_energy_site(int index) {
    
    //get position in lattice 
    int ix = index%nx;
    int iy = index/nx;
    int n_cell_x = nx;
    int n_cell_y = nx;
    
    // get neighbours with pbc
    int nb[4], e;
    nb[0] = ((ix + 1) % n_cell_x + ((iy + n_cell_y) % n_cell_y) * n_cell_x); //B
    nb[1] = ((ix + n_cell_x) % n_cell_x + ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x); //F
    nb[2] = ((ix + n_cell_x - 1) % n_cell_x + ((iy + n_cell_y) % n_cell_y) * n_cell_x); //H
    nb[3] = ((ix + n_cell_x) % n_cell_x + ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x); //D
    
    //calc energy
    e = -1*lattice[index]*(lattice[nb[0]] + lattice[nb[1]] + lattice[nb[2]] + lattice[nb[3]]);
    return e;
}

int get_total_energy() {
    int tote=0;
    for(int i=0; i<N; i++) tote += get_energy_site(i);
    return tote;
}

int get_total_magnetization() {
    int totm=0;
    for(int i=0; i<N; i++) totm += lattice.at(i);
    return totm;
}

//*********** trial move ***********//
bool test_flip(int index, int &de) {
    de = -2*get_energy_site(index);

    if(de < 0) 
        return true;
    
    else if(dist(rng) < std::exp(-de/T)) 
        return true;
    
    else 
        return false;
}

bool test_flip(int de) {
    if(de <= 0) 
        return true;
    else if(dist(rng) < std::exp(-de/T)) 
        return true;
    else 
        return false;
}

// transient steps [Metropolis ONLY]
void discard_before_eq(int mcs) {
    int de;
    
    //Monte Carlo loop
    for(long i=0; i<mcs; i++) {
        
        //Metropolis loop
        for(int i=0; i<N; i++) {
            int rand_pos = static_cast<int>(dist(rng)*N);
            
            if(test_flip(rand_pos,de)) {
                //flip 
                lattice[rand_pos] = -lattice[rand_pos];
            }
            
        } // Metropolis loop end
        
    } //Monte Carlo loop end 
}

//*********** main ***********//
int main(int argc, char *argv[]) {
    //set parameters
    nx = atoi(argv[1]);
    N = nx*nx;
    T = atof(argv[2]);
    mcs = atoi(argv[3]);
    int bin_count = atoi(argv[4]);
    fthermo.open(argv[5], std::ios::out);
    string fname =(string(argv[5])+".debug");
    fthermo2.open(fname.c_str(), std::ios::out);
    
    //thermos
    double E,M,m;
    double etot,mtot;
    
    cout << "lattice -> " << nx << "x" << nx << ", T=" << T << ", mcs=" << mcs << endl;
    cout << "nbins -> " << bin_count << "thermo_file -> " << argv[5] << endl;
    
    //SUS bins
    //vector<long> sus_hist;
    vector<double> P;
    vector<int> lattice_saved;
    
    //int window_count = bin_count-1;///2;
    //double bin_min = -N;
    //double bin_max = N;
    //double bin_width = (bin_max-bin_min)/double(bin_count);
    cout << "N -> " << N << endl;
    
    //sus_hist.resize(bin_count);
    P.resize(N+1);
    
    
    //for(auto &i:sus_hist) i=0;
    
    //open files
    flat.open("lattice.dat", std::ios::out);
    
    //allocate
    lattice.resize(N);
    lattice_saved.resize(N);
    
    //Generate
    generate_lattice(-1);
    M = get_total_magnetization();
    
    //print_lattice();
    //cout << get_total_energy() << " " << get_total_magnetization() << endl;
    
    //DEBUG BINS
    //for(int i=0; i<bin_count; i++) {
    //    int bin_mag_min = i    *bin_width+bin_min;
    //    int bin_mag_max = (i+1)*bin_width+bin_min;
    //    cout << i << ": " << bin_mag_min << " -> " << bin_mag_max << endl;
    //}
    cout << "------------------------------------" << endl;
    
    int ibin=0;
    int window_edges[2];
    int l,r,de;
    P[0]=1;
    double Nd=0;
    lattice_saved = lattice;
    int lattice_saved_count = 0;
    int M_trial;
    int rand_pos;
    
    for(int i=0; i<N; i++) { //dir=0 => left; dir=1 =>right
        cout << i << "\r";
        cout.flush();
        
        lattice = lattice_saved;
        //lattice.swap(lattice_saved);
        
        l=0;
        r=0;
        lattice_saved_count = 0;
        
        M = get_total_magnetization();
        
        // Monte Carlo loop
        for(long j=0; j<mcs; j++) {
            
            // calculate change in energy & magnetization 
            // if a spin at random site rand_pos is flipped
            
            rand_pos = int(dist(rng)*N); 
            M_trial = M - 2*lattice[rand_pos] ;
                        
            //***************** SUS ******************
            int pre  = -N+2*i;
            int post = -N+2*(i+1);
            
            // IF M_trial is in the current window
            if(M_trial >= pre && M_trial <= post ) {
                
                //now lets see whether this move is ACCEPTED by Metropolis                    
                if(test_flip(rand_pos,de)) {
                    lattice[rand_pos] *= -1; //FLIP
                    M = M_trial; // As the move is accepted, so our new M is equal to M_trial
                } 
                
                if(M == pre) {
                    l++;
                }
                
                if(M == post) {
                    r++;
                    //if(lattice_saved_count==0) {
                        
                        // ?? What is the FASTEST way to copy a c++ vector to another one ??
                        //lattice_saved = lattice;
                        lattice_saved.assign(lattice.begin(), lattice.end());
                        lattice_saved_count++;
                    //}
                }
            } 
            
            //If this move is outside the window => REJECT it
            else {
                if(M_trial < pre)  l++;
                if(M_trial > post) r++;                
            }
            
        } //Monte Carlo loop end
        
        
        P[i + 1] = P[i] + log(r*pow(l,-1));
        Nd += exp(P[(i + 1)]);
        cout << l << " " << r << " " << Nd << " " << lattice_saved_count << endl;
    
    } // SUS loop ends
    
    
    // op to file
    for (int i = 0; i < N; i++) {
        fthermo << -N+2*i << " " << exp(P[i + 1])/Nd << endl;
    }
    
    // close all
    flat.close();
    fthermo.close();
    
    return 0;
}



// 7207101844
