//
// Created by Vitto Resnick on 3/29/24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

// Check if a number is an integer
bool is_integer(double k){
    return floor(k) == k;
}

// Factorial function
int factorial(int n){
    int res = 1,i;
    for (i=1;i<=n;i++){
        res *= i;
    }
    return res;
}

// Double factorial
int dfact(int n){
    int i;double res=1.0;
    for(i=n;i>=1;i-=2){
        res *=i;
    }
    return res;
}

// Binomial Coefficient Example: m choose n
int binomCoeff(int m, int n){
    return factorial(m)/(factorial(n)*factorial(m-n));
}

// Calculation one of three directional components for SAB (overlap integral of prod of 2 gauss)
double SxABComp(double XA, double XB, double alpha, double beta, double lA, double lB){

    double P  = exp(-alpha*beta*pow(XA-XB,2)/(alpha + beta)); // Calculate prefactor
    double XP = (alpha*XA + beta*XB)/(alpha + beta);
    double doubleSum = 0; // Initialize double sum

    // Compute  double sound
    for     (int i = 0; i < lA+1; i++){
        double innerSum = 0;
        for (int j = 0; j < lB+1; j++){
            if ((i+j)% 2 == 0){ // Only do even i+j terms
                double summand = binomCoeff(lA,i)*binomCoeff(lB,j)*dfact(i+j-1)*pow(XP-XA,lA-i)*pow(XP-XB,lB-j)/pow(2*(alpha+beta),(i+j)/2);
                innerSum += summand;
            }
        }
        doubleSum += innerSum;
    }
    return P*sqrt(M_PI/(alpha + beta))*doubleSum;
}

// Find normalization constant for a given k primitive gauss
double NormConst(double X, double Y, double Z, double alpha_k, double l, double m, double n){
    double SAA = 1; // Initialize SAA product
    SAA *= SxABComp(X,X,alpha_k,alpha_k,l,l); // Compute SxAA
    SAA *= SxABComp(Y,Y,alpha_k,alpha_k,m,m); // Compute SyAA
    SAA *= SxABComp(Z,Z,alpha_k,alpha_k,n,n); // Compute SzAA

    // normalization constants are defined such that the overlap of each primitive with itself is equal to one
    double N_k_lmn = 1/sqrt(SAA);
    return N_k_lmn;
}

// Process basis data from basis files
vector<vec> ProcessBasisData(string H_STO3G, string C_STO3G, string N_STO3G, string O_STO3G, string F_STO3G){
    double H_alpha_1s_1,H_alpha_1s_2,H_alpha_1s_3;
    double H_d_1s_1,H_d_1s_2,H_d_1s_3;

    double C_alpha_2s_2p_1,C_alpha_2s_2p_2,C_alpha_2s_2p_3;
    double C_d_2s_1       ,C_d_2s_2       ,C_d_2s_3       ;
    double C_d_2p_1       ,C_d_2p_2       ,C_d_2p_3       ;

    double N_alpha_2s_2p_1,N_alpha_2s_2p_2,N_alpha_2s_2p_3;
    double N_d_2s_1       ,N_d_2s_2       ,N_d_2s_3       ;
    double N_d_2p_1       ,N_d_2p_2       ,N_d_2p_3       ;

    double O_alpha_2s_2p_1,O_alpha_2s_2p_2,O_alpha_2s_2p_3;
    double O_d_2s_1       ,O_d_2s_2       ,O_d_2s_3       ;
    double O_d_2p_1       ,O_d_2p_2       ,O_d_2p_3       ;

    double F_alpha_2s_2p_1,F_alpha_2s_2p_2,F_alpha_2s_2p_3;
    double F_d_2s_1       ,F_d_2s_2       ,F_d_2s_3       ;
    double F_d_2p_1       ,F_d_2p_2       ,F_d_2p_3       ;


    ifstream inputFile0(H_STO3G);
    inputFile0 >> H_alpha_1s_1     >> H_d_1s_1;
    inputFile0 >> H_alpha_1s_2     >> H_d_1s_2;
    inputFile0 >> H_alpha_1s_3     >> H_d_1s_3;
    inputFile0.close();

    ifstream inputFile1(C_STO3G);
    inputFile1 >> C_alpha_2s_2p_1 >> C_d_2s_1 >> C_d_2p_1;
    inputFile1 >> C_alpha_2s_2p_2 >> C_d_2s_2 >> C_d_2p_2;
    inputFile1 >> C_alpha_2s_2p_3 >> C_d_2s_3 >> C_d_2p_3;
    inputFile1.close();

    ifstream inputFile2(N_STO3G);
    inputFile2 >> N_alpha_2s_2p_1 >> N_d_2s_1 >> N_d_2p_1;
    inputFile2 >> N_alpha_2s_2p_2 >> N_d_2s_2 >> N_d_2p_2;
    inputFile2 >> N_alpha_2s_2p_3 >> N_d_2s_3 >> N_d_2p_3;
    inputFile2.close();

    ifstream inputFile3(O_STO3G);
    inputFile3 >> O_alpha_2s_2p_1 >> O_d_2s_1 >> O_d_2p_1;
    inputFile3 >> O_alpha_2s_2p_2 >> O_d_2s_2 >> O_d_2p_2;
    inputFile3 >> O_alpha_2s_2p_3 >> O_d_2s_3 >> O_d_2p_3;
    inputFile3.close();

    ifstream inputFile4(F_STO3G);
    inputFile4 >> F_alpha_2s_2p_1 >> F_d_2s_1 >> F_d_2p_1;
    inputFile4 >> F_alpha_2s_2p_2 >> F_d_2s_2 >> F_d_2p_2;
    inputFile4 >> F_alpha_2s_2p_3 >> F_d_2s_3 >> F_d_2p_3;
    inputFile4.close();
    vector<vec> data;
    data.push_back(vec{H_alpha_1s_1,H_alpha_1s_2,H_alpha_1s_3});
    data.push_back(vec{H_d_1s_1,H_d_1s_2,H_d_1s_3});

    data.push_back(vec{C_alpha_2s_2p_1,C_alpha_2s_2p_2,C_alpha_2s_2p_3});
    data.push_back(vec{C_d_2s_1       ,C_d_2s_2       ,C_d_2s_3       });
    data.push_back(vec{C_d_2p_1       ,C_d_2p_2       ,C_d_2p_3       });

    data.push_back(vec{N_alpha_2s_2p_1,N_alpha_2s_2p_2,N_alpha_2s_2p_3});
    data.push_back(vec{N_d_2s_1       ,N_d_2s_2       ,N_d_2s_3       });
    data.push_back(vec{N_d_2p_1       ,N_d_2p_2       ,N_d_2p_3       });

    data.push_back(vec{O_alpha_2s_2p_1,O_alpha_2s_2p_2,O_alpha_2s_2p_3});
    data.push_back(vec{O_d_2s_1       ,O_d_2s_2       ,O_d_2s_3       });
    data.push_back(vec{O_d_2p_1       ,O_d_2p_2       ,O_d_2p_3       });

    data.push_back(vec{F_alpha_2s_2p_1,F_alpha_2s_2p_2,F_alpha_2s_2p_3});
    data.push_back(vec{F_d_2s_1       ,F_d_2s_2       ,F_d_2s_3       });
    data.push_back(vec{F_d_2p_1       ,F_d_2p_2       ,F_d_2p_3       });
    return data;
}

// Select exp, contraction, quantum num, and coord constants for a given basis function
tuple<vec,vec,vec,vec,vec> FindConsts(vector<vec> basis_data, int atom, vec R_center,string orbital){
    // Exponent Data
    vec H_alpha_1s    = basis_data[0];
    vec X_alpha_2s_2p;
    // Contraction Coefficient Data
    vec H_d_1s = basis_data[1];
    vec X_d_2s;
    vec X_d_2p;
    if (atom == 6){
        X_alpha_2s_2p = basis_data[2];
        X_d_2s = basis_data[3];
        X_d_2p = basis_data[4];
    } else if (atom == 7){
        X_alpha_2s_2p = basis_data[5];
        X_d_2s = basis_data[6];
        X_d_2p = basis_data[7];
    } else if (atom == 8){
        X_alpha_2s_2p = basis_data[8];
        X_d_2s = basis_data[9];
        X_d_2p = basis_data[10];
    } else if (atom == 9){
        X_alpha_2s_2p = basis_data[11];
        X_d_2s = basis_data[12];
        X_d_2p = basis_data[13];
    }

    // Quantum Numbers
    const vec lms_s = vec({0,0,0});
    const vec lms_px = vec({1,0,0});
    const vec lms_py = vec({0,1,0});
    const vec lms_pz = vec({0,0,1});

    // Initialize Output Vectors
    vec exponents;
    vec contraCoeffs;
    vec quantNums;
    if (atom == 1){
        exponents    = H_alpha_1s;
        contraCoeffs = H_d_1s;
        quantNums    = lms_s;
    } else if (atom == 6 || atom == 7 || atom == 8 || atom == 9){
        if        (orbital == "2s" ){
            exponents    = X_alpha_2s_2p ;
            contraCoeffs = X_d_2s;
            quantNums    = lms_s;
        } else if (orbital == "2px"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_px;
        } else if (orbital == "2py"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_py;
        } else if (orbital == "2pz"){
            exponents    = X_alpha_2s_2p;
            contraCoeffs = X_d_2p;
            quantNums    = lms_pz;
        }
    }
    // For each basis function, ωμ(r) (for μ = 1 · · · N ) there will be ...

    // (i) a center, R,
    double X = R_center(0);
    double Y = R_center(1);
    double Z = R_center(2);

    // (ii) 3 quantum numbers, (l , m, n),
    double l = quantNums(0);
    double m = quantNums(1);
    double n = quantNums(2);

    // and (iii) information about 3 primitive functions:
    // 3 exponents, αk,
    // 3 corresponding contraction coefficients, dk, and
    // 3 normalization constants, N_k_lmn
    int K = contraCoeffs.size();
    vec normConsts(K, fill::ones);
    for (int k = 0; k < K; k++){
        double alpha_k = exponents(k);
        double N_k_lmn = NormConst(X,Y,Z,alpha_k,l,m,n);
        normConsts(k) = N_k_lmn;
    }

    return make_tuple(R_center,quantNums,exponents,contraCoeffs,normConsts);
}

// Read in input file
tuple<int,vec,int,int,vector<vector<double>>,vector<int>,vector<vector<double>>,vector<int>,vector<int>> read_input_file(string file_name){
    // Read in the coordinates, in the format: E X Y Z for each atom, where E is the element (handle at least H and C).

    // Read in input file
    ifstream inputFile(file_name);

    // Initialize vars
    int n_atoms;          // Initialize total number of atoms
    int charge;
    inputFile >> n_atoms >> charge; // Set total number of atoms and charge
    const double Bohr_A = 0.52917706; // angstroms in 1 bohr
    int a = 0; // Number carbons
    int b = 0;
    vector<vector<double>> xyz_list;       // Initialize list for atoms' xyz coordinates
    vector<int> atom_list;                 // Initialize list for atoms' identities
    vector<vector<double>> basis_xyz_list; // Initialize list for atoms' xyz coordinates
    vector<int> basis_atom_list;           // Initialize list for atoms' identities
    vector<int> basis_atom_nums;
    vec Z = vec(n_atoms,fill::zeros);
    // Read in atom identity and xyz coordinates
    for (int i = 0; i < n_atoms; ++i) {  // Iterate through every atom
        int atom;
        double x, y, z;                    // Initialize atom identity and xyz coordinates
        inputFile >> atom >> x >> y >> z ; // Set atomic number/atom identity and xyz coordinates
        //x = x/Bohr_A;
        //y = y/Bohr_A;
        //z = z/Bohr_A;
        if (atom != 1 && atom != 6 && atom != 7 && atom != 8 && atom != 9) {      // If a given atom is not H or C, throw an error
            cerr << "Atom No." << i+1 << ": This atom is not a hydrogen, carbon, nitrogen, oxygen, fluorine!" << endl;
        } else if (atom == 6 || atom == 7 || atom == 8 || atom == 9){
            a = a + 1;
            basis_atom_list.insert(basis_atom_list.end(), {atom,atom,atom,atom});
            basis_xyz_list.insert(basis_xyz_list.end(), { {x, y, z},{x, y, z},{x, y, z},{x, y, z} });
            basis_atom_nums.insert(basis_atom_nums.end(),{i,i,i,i});
        } else if (atom == 1){
            b = b + 1;
            basis_atom_list.push_back(atom);
            basis_xyz_list.push_back({x, y, z});
            basis_atom_nums.push_back(i);
        }

        atom_list.push_back(atom);         // Append this atom's atomic number/atom identity to list
        xyz_list.push_back({x, y, z});     // Append this atom's xyz coordinates to list

        // ZA, valence atomic number of
        if (atom == 1){
            Z(i)=1;
        } else if (atom == 6){
            Z(i)=4;
        } else if (atom == 7){
            Z(i)=5;
        } else if (atom == 8){
            Z(i)=6;
        } else if (atom == 9){
            Z(i)=7;
        }
    }
    inputFile.close();                     // Close the txt file

    tuple<int,vec,int,int,vector<vector<double>>,vector<int>,vector<vector<double>>,vector<int>,vector<int>> output = make_tuple(n_atoms,Z,a,b,xyz_list,atom_list,basis_xyz_list,basis_atom_list,basis_atom_nums);
    return output;
}

// Build Basis Functions
tuple<vector<tuple<vec,vec,vec,vec,vec>>,vector<string>> build_basis(vector<vec> basis_data, int N, vector<int> basis_atom_list, vector<vector<double>> basis_xyz_list){
    //input basis_data, N, basis_atom_list, basis_xyz_list
    // Build a list of the basis functions, which are contracted gaussians,
    // Iterate through N AOs and construct basis functions
    vector<string> row2_orbital_bank = {"2s","2px","2py","2pz"};
    vector<tuple<vec,vec,vec,vec,vec>> basis_func_constants;
    vector<string> basis_orbital_list;

    for (int i = 0; i < N; i++){

        int atom = basis_atom_list[i]; // select atom
        vector<double> center = basis_xyz_list[i]; // find coordinates

        string orbital;
        if (atom == 1){
            orbital = "1s"; // If it's H, only 1s orbital
        } else if (atom == 6 || atom == 7 || atom == 8 || atom == 9){ // If it's C, for each C, there's 2s, 2px, 2py, and 2pz
            orbital = row2_orbital_bank[0];
            row2_orbital_bank.erase(row2_orbital_bank.begin());
        }
        if (row2_orbital_bank.size() == 0){
            row2_orbital_bank = {"2s","2px","2py","2pz"};
        }

        basis_orbital_list.push_back(orbital); // Collect basis functions' orbitals + find their characteristic consts
        basis_func_constants.push_back(FindConsts(basis_data,atom, center, orbital));

    }
    return make_tuple(basis_func_constants,basis_orbital_list);
}

// S, AO basis overlap integral matrix
mat build_s(int N, vector<tuple<vec,vec,vec,vec,vec>> basis_func_constants){
    // Find all normalization constants for the basis functions
    mat All_Normalization_Constants(3,N,fill::zeros);

    // Run a loop over all your basis functions
    for (int i = 0; i < N; i++){
        // get the normalization constants for the 3 primitives that make up each basis function
        auto[R_center,quantNums,exponents,contraCoeffs,normConsts] = basis_func_constants[i];
        // and save them in an array.
        All_Normalization_Constants.col(i) = normConsts;
    }


    // Contracted overlap integral S
    mat S(N,N,fill::zeros);

    // Iterate over pairs of basis functions
    for (int mu = 0; mu < N; mu++){
        for (int nu = 0; nu < N; nu++){
            // For this given pair of basis functions
            // R_center,quantNums,exponents,contraCoeffs,normConsts for each atom in the pair
            auto[R_mu,lmn_mu,a_mu,d_mu,N_mu] = basis_func_constants[mu];
            auto[R_nu,lmn_nu,a_nu,d_nu,N_nu] = basis_func_constants[nu];

            // # of contracted gaussians
            int K = d_mu.size();
            int L = d_nu.size();

            double S_mu_nu = 0;

            for (int k = 0; k < K; k++){
                for (int l = 0; l < L; l++){
                    double d_k_mu     = d_mu(k);
                    double d_l_nu     = d_nu(l);
                    double N_k_mu     = N_mu(k);
                    double N_l_nu     = N_nu(l);
                    double alpha_k_mu = a_mu(k);
                    double alpha_l_nu = a_nu(l);

                    double Xm = R_mu(0);
                    double Ym = R_mu(1);
                    double Zm = R_mu(2);

                    double Xn = R_nu(0);
                    double Yn = R_nu(1);
                    double Zn = R_nu(2);

                    double lm = lmn_mu(0);
                    double mm = lmn_mu(1);
                    double nm = lmn_mu(2);

                    double ln = lmn_nu(0);
                    double mn = lmn_nu(1);
                    double nn = lmn_nu(2);

                    double S_k_l = 1; // Initialize SAB product
                    S_k_l *= SxABComp(Xm,Xn,alpha_k_mu,alpha_l_nu,lm,ln); // Compute Sxkl
                    S_k_l *= SxABComp(Ym,Yn,alpha_k_mu,alpha_l_nu,mm,mn); // Compute Sykl
                    S_k_l *= SxABComp(Zm,Zn,alpha_k_mu,alpha_l_nu,nm,nn); // Compute Szkl

                    S_mu_nu += d_k_mu * d_l_nu * N_k_mu * N_l_nu * S_k_l;
                }
            }
            S(mu,nu) = S_mu_nu;
        }
    }
    return S;
}

tuple<double,double,double,double,double,double,double,double>gamma_atom_constants(vector<vec> basis_data,vector<int> atom_list, vec RA, vec RB, int A, int B, int k, int kp, int l, int lp){
    vec d_A;
    vec d_B;
    vec alpha_A;
    vec alpha_B;

    int atom = atom_list[A];
    vec d;
    vec alpha;
    if (atom == 1){
        alpha = basis_data[0];
        d = basis_data[1];
    } else if (atom == 6){
        alpha = basis_data[2];
        d = basis_data[3];
    } else if (atom == 7){
        alpha = basis_data[5];
        d = basis_data[6];
    } else if (atom == 8){
        alpha = basis_data[8];
        d = basis_data[9];
    } else if (atom == 9){
        alpha = basis_data[11];
        d = basis_data[12];
    }
    d_A = d;
    alpha_A = alpha;

    atom = atom_list[B];
    if (atom == 1){
        alpha = basis_data[0];
        d = basis_data[1];
    } else if (atom == 6){
        alpha = basis_data[2];
        d = basis_data[3];
    } else if (atom == 7){
        alpha = basis_data[5];
        d = basis_data[6];
    } else if (atom == 8){
        alpha = basis_data[8];
        d = basis_data[9];
    } else if (atom == 9){
        alpha = basis_data[11];
        d = basis_data[12];
    }
    d_B = d;
    alpha_B = alpha;

    // Extract contraction coefficients
    double d_k_sA  = d_A(k);
    double d_kp_sA = d_A(kp);
    double d_l_sB  = d_B(l);
    double d_lp_sB = d_B(lp);

    // Extract exponents
    double alpha_k  = alpha_A(k);
    double alpha_kp = alpha_A(kp);
    double alpha_l  = alpha_B(l);
    double alpha_lp = alpha_B(lp);

    // Calculate Normalization Constants
    double N_k_sA  = NormConst(RA(0), RA(1), RA(2), alpha_k , 0, 0, 0);
    double N_kp_sA = NormConst(RA(0), RA(1), RA(2), alpha_kp, 0, 0, 0);
    double N_l_sB  = NormConst(RB(0), RB(1), RB(2), alpha_l , 0, 0, 0);
    double N_lp_sB = NormConst(RB(0), RB(1), RB(2), alpha_lp, 0, 0, 0);

    // Calculate d primes
    double dp_k_sA  = d_k_sA *N_k_sA;
    double dp_kp_sA = d_kp_sA*N_kp_sA;
    double dp_l_sB  = d_l_sB *N_l_sB;
    double dp_lp_sB = d_lp_sB*N_lp_sB;

    return make_tuple(alpha_k,alpha_kp,alpha_l,alpha_lp,dp_k_sA,dp_kp_sA,dp_l_sB,dp_lp_sB);
}

// Build Gamma = all two-center 2-electron repulsion integrals evaluated over the square of the valence s
// orbital centered on atoms A and B
mat build_g(int n_atoms, vector<vector<double>> xyz_list,vector<vec> basis_data,vector<int> atom_list){
    // Contracted overlap integral S

    mat gammaAB(n_atoms,n_atoms,fill::zeros);

    for (int A = 0; A < n_atoms; A++){
        for (int B = 0; B < n_atoms; B++){

            //Define RA and RB
            vector<double> RA_raw = xyz_list[A];
            vector<double> RB_raw = xyz_list[B];
            vec RA = vec{RA_raw[0],RA_raw[1],RA_raw[2]};
            vec RB = vec{RB_raw[0],RB_raw[1],RB_raw[2]};

            double g = 0;

            // sum over two-electron integrals over primitive Gaussians
            for             (int k = 0; k < 3; k++){
                for         (int kp = 0; kp < 3; kp++){
                    for     (int l = 0; l < 3; l++){
                        for (int lp = 0; lp < 3; lp++){
                            double summand;

                            auto[alpha_k,alpha_kp,alpha_l,alpha_lp,dp_k_sA,dp_kp_sA,dp_l_sB,dp_lp_sB]= gamma_atom_constants(basis_data,atom_list,RA,RB,A,B,k,kp,l,lp);

                            // Compute 0^0 integral
                            // 6-dimensional 2-center integral over the product of 4 primitive s-type Gaussian functions
                            double sigma_A = 1/(alpha_k+alpha_kp);
                            double sigma_B = 1/(alpha_l+alpha_lp);

                            double UA = pow((M_PI*sigma_A),1.5);
                            double UB = pow((M_PI*sigma_B),1.5);

                            double V2 = 1/(sigma_A+sigma_B);

                            double diff2 = pow(norm(RA-RB),2);
                            double T = V2*diff2;

                            double oo = 0;
                            if (T==0){
                                oo = UA*UB*sqrt(2*V2)*sqrt(2.0/M_PI);
                            } else {
                                oo = UA*UB*sqrt(1.0/diff2)*erf(sqrt(T));
                            }

                            summand = dp_k_sA*dp_kp_sA*dp_l_sB*dp_lp_sB*oo;
                            g += summand;
                        }
                    }
                }
            }
            gammaAB(A,B) = 27.2114079527*g; //convert to eV from atomic units
        }
    }
    return gammaAB;
}

// Build core hamiltonian matrix, which is independent of the density matrix
mat build_h(int N, int n_atoms, vector<int> basis_atom_nums, vector<int> basis_atom_list, vector<string> basis_orbital_list, vec Z, mat IAb, mat gammaAB, mat s){
    mat H(N,N,fill::zeros);
    // inputs basis_atom_nums, basis_orbital_list, IAb, n_atoms, p_tot, Z, gamma, p, s
    // Iterate over pairs of basis functions

    for (int mu = 0; mu < N; mu++) {
        for (int nu = 0; nu < N; nu++) {
            double h;
            if (mu == nu){
                // calculate summation
                int A = basis_atom_nums[mu];

                double IA;
                string orbital = basis_orbital_list[mu];
                if (orbital == "1s" || orbital == "2s"){
                    IA = IAb(0,basis_atom_list[mu]-1);
                } else if (orbital == "2px" || orbital == "2py" || orbital == "2pz"){
                    IA = IAb(1,basis_atom_list[mu]-1);
                }

                double summation = 0;
                for (int B = 0; B < n_atoms; B++){
                    if (B != A){
                        summation += Z(B)*gammaAB(A,B);
                    }
                }
                h = -IA-(Z(A)-0.5)*gammaAB(A,A)-summation;
            } else {
                int A = basis_atom_nums[mu];
                int B = basis_atom_nums[nu];

                double nbetaA = IAb(2,basis_atom_list[mu]-1);
                double nbetaB = IAb(2,basis_atom_list[nu]-1);

                h = -0.5*(nbetaA+nbetaB)*s(mu,nu);
            }
            H(mu,nu) = h;
        }
    }
    return H;
}

// Build density matrix for p alpha electrons or q beta electrons
mat build_p(mat C_a, int p) {

    int n_rows = C_a.n_rows; // since the size of a matrix is the total # of elements, and c_alpha is always square

    // Initialize appropriately sized matrix
    mat C_occupied(n_rows, p, arma::fill::zeros);

    // slice the c_alpha square matrix into a rectangular matrix, just selecting the first p columns
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < p; j++) {
            C_occupied(i,j) = C_a(i,j);
        }
    }
    // occupied coefficient matrix times transpose
    mat P_mat = C_occupied*C_occupied.t();
    return P_mat;
}

// Total density on each atom
vec build_pAA(int N, int n_atoms, vector<int> basis_atom_nums, mat p_alpha, mat p_beta){
    mat p_tot = p_alpha + p_beta;
    vec pAA(n_atoms,fill::zeros);
    for     (int A = 0; A < n_atoms; A++){
        double pAA_A = 0;
        for (int mu = 0; mu < N; mu++) {
            if (basis_atom_nums[mu]==A){
                pAA_A += p_tot(mu,mu);
            }
        }
        pAA(A) = pAA_A;
    }
    return pAA;
}

// Assemble the matrix elements of the CNDO/2 Fock operator, f (in the AO basis)
mat build_f(int N, int n_atoms, vector<int> basis_atom_nums, vector<int> basis_atom_list, vector<string> basis_orbital_list, vec Z, mat IAb, mat gammaAB, mat s, mat p, vec p_tot){

    mat F(N,N,fill::zeros);
    // inputs basis_atom_nums, basis_orbital_list, IAb, n_atoms, p_tot, Z, gamma, p, s
    // Iterate over pairs of basis functions
    for (int mu = 0; mu < N; mu++) {
        for (int nu = 0; nu < N; nu++) {
            double f;

            if (mu == nu){
                // Diagonal Fock matrix elements
                // Do unrestricted equations

                // calculate summation
                int A = basis_atom_nums[mu];

                // Semi-empirical parameters of ionization energies and electron affinities
                double IA;
                string orbital = basis_orbital_list[mu];
                if (orbital == "1s" || orbital == "2s"){
                    IA = IAb(0,basis_atom_list[mu]-1);
                } else if (orbital == "2px" || orbital == "2py" || orbital == "2pz"){
                    IA = IAb(1,basis_atom_list[mu]-1);
                }

                double summation = 0;
                for (int B = 0; B < n_atoms; B++){
                    if (B != A){
                        summation += (p_tot(B)-Z(B)) *gammaAB(A,B);
                    }
                }
                f = -IA+((p_tot(A)-Z(A))-(p(mu,mu)-0.5))*gammaAB(A,A)+summation;
            } else {
                // Off-diagonal Fock matrix elements
                // Do unrestricted equations

                int A = basis_atom_nums[mu];
                int B = basis_atom_nums[nu];

                double nbetaA = IAb(2,basis_atom_list[mu]-1);
                double nbetaB = IAb(2,basis_atom_list[nu]-1);
                f = -0.5*(nbetaA + nbetaB) * s(mu,nu) - p(mu,nu) * gammaAB(A,B);
            }
            F(mu,nu) = f;
        }
    }
    return F;
}

// Fock Matrix to MO Coefficients
mat F_to_C(mat F) {
    // Find MO Coeffs
    vec energies;
    mat C;
    eig_sym(energies, C, F);

    return C;
}

// Fock Matrix to Energy eigen values
mat F_to_ep(mat F) {
    vec energies;
    mat C;
    eig_sym(energies, C, F); // can use eig_sym because F is guaranteed to be symmetric
    return energies;
}

// Electron energy
double E_electron(int N, mat H, mat P_a, mat P_b, mat F_a, mat F_b){
    double E = 0;
    for (int mu = 0; mu < N; mu++){
        for (int nu = 0; nu < N; nu++){
            E += P_a(mu,nu)*(H(mu,nu)+F_a(mu,nu))+P_b(mu,nu)*(H(mu,nu)+F_b(mu,nu));
        }
    }
    return E/2; //avoid double counting
}

// Nuclear Repulsion Energy
double E_nuc_rep(int n_atoms,vec Z, vector<vector<double>> xyz_list){
    double E = 0;
    for (int A = 0; A < n_atoms; A++){
        for (int B = 0; B < n_atoms; B++){
            if (B > A){ // Don't count the same atom
                vector<double> RA_raw = xyz_list[A];
                vector<double> RB_raw = xyz_list[B];
                double xA = RA_raw[0], yA = RA_raw[1], zA = RA_raw[2];
                double xB = RB_raw[0], yB = RB_raw[1], zB = RB_raw[2];

                double R_AB = sqrt(pow(xA-xB,2) + pow(yA-yB,2) + pow(zA-zB,2));

                E += Z(A)*Z(B)/R_AB;
            }
        }
    }
    return E*27.2114079527; //convert to eV and avoid double counting
}

// Vector norm calculator function
double Rij(double xi, double xj, double yi, double yj, double zi, double zj) {
    return sqrt(pow(xi - xj, 2) + pow(yi - yj, 2) + pow(zi - zj, 2));
}

// Convert 2D vector vector<vector<double>> to matrix
mat VecVec2Mat(int n, vector<vector<double>> xyz_list){
    mat new_matrix(3,n,fill::zeros);
    for (int i = 0; i < n; ++i){
        for (int j=0; j < 3; ++j){
            new_matrix(j,i)=xyz_list[i][j];
        }
    }
    return new_matrix;
}

// Convert matrix to 2D vector vector<vector<double>>
vector<vector<double>> Mat2VecVec(int n, mat coord_mat){
    vector<vector<double>> new_vec;
    for (int i = 0; i < n; ++i){
        vector<double> coord = {coord_mat(0,i),coord_mat(1,i),coord_mat(2,i)};
        new_vec.push_back(coord);
    }
    return new_vec;
}

// Vectorize Gradient -> Calculate Norm value of Gradient
double GradNorm(int n,mat Gradient){
    return norm(vectorise(Gradient));
}

// Calculate Unit Vector (matrix) of Gradient or Normalize Gradient
mat UnitGrad(int n,mat Gradient){
    return Gradient/GradNorm(n,Gradient);
}

tuple<mat,mat,vec,mat,mat,mat,mat> CNDO_converge(int N, int n_atoms, vector<int> basis_atom_nums, vector<int> basis_atom_list, mat S, mat gammaAB, mat IAb, vector<string> basis_orbital_list, vec Z, int p, int q){

    int iterator = 0;
    bool convergence = false;
    double tol = 1e-8;
    mat P_a, P_b, P_a_old, P_b_old;
    mat C_a, C_b, C    ;
    mat F_a, F_b       ;
    vec p_tot;


    while (!convergence){
        //cout << "Iteration: " << iterator << endl;

        // Given the AO basis total density matrix, p, where,
        // for an open shell molecule with p α electrons and q β electrons,
        if (iterator == 0){
            // Base Case
            // 1. Guess pα = pβ = 0, Take an initial guess of the density matrix
            P_a = mat(N,N,fill::zeros);
            P_b = mat(N,N,fill::zeros);
            p_tot = build_pAA(N,n_atoms,basis_atom_nums,P_a,P_b); // Total densities on each atom
        }

        // Assemble the matrix elements of the CNDO/2 Fock operator, f (in the AO basis)
        // 2. Build Fock matrix fα and fβ.
        F_a = build_f(N,n_atoms,basis_atom_nums,basis_atom_list,basis_orbital_list,Z,IAb,gammaAB,S,P_a,p_tot);
        F_b = build_f(N,n_atoms,basis_atom_nums,basis_atom_list,basis_orbital_list,Z,IAb,gammaAB,S,P_b,p_tot);

        //F_a.print("Fa");
        //F_b.print("Fb");

        // 3. Solve the eigenvalue problems, Eqs. 2.3 and 2.4, to obtain new MO coefficients, cα and cβ
        C_a = F_to_C(F_a);
        C_b = F_to_C(F_b);

        //cout << "after solving eigenvalue equation:" << iterator << endl;

        //C_a.print("Ca");
        //C_b.print("Cb");

        // 4. Copy the old density matrices to pα old and pβ old
        P_a_old = P_a;
        P_b_old = P_b;

        // 5. Assemble new density matrices by occupying the p lowest energy α MOs and the q
        // lowest energy β MOs and using Eqs. 1.1, 1.2 and 1.3
        //cout << "p = " << p << " q = " << q << endl;
        P_a = build_p(C_a,p);
        P_b = build_p(C_b,q);

        //P_a.print("Pa_new");
        //P_b.print("Pb_new");
        p_tot = build_pAA(N,n_atoms,basis_atom_nums,P_a,P_b); // Total densities on each atom
        //p_tot.print("P_t");

        iterator += 1;
        // 6. If the maximum magnitude of the change in the α and β density matrices is less than
        // a tolerance you specify (e.g. 10−6), then you have converged. Otherwise return to step 2.
        convergence = approx_equal(P_a, P_a_old, "absdiff", tol);
    }
    return make_tuple(P_a,P_b,p_tot,F_a,F_b,C_a,C_b);
}

// Evaluate the CNDO/2 energy
double CNDO_energy(mat F_a, mat F_b, int N, int n_atoms, mat H, mat P_a, mat P_b, vec Z, vector<vector<double>> xyz_list, string file_name, bool print_energy){
    //cout << "\nFinal Converged Molecule:" << endl;
    // 3. corresponding eigenvalues, εα and εβ
    vec E_a = F_to_ep(F_a);
    vec E_b = F_to_ep(F_a);

    //E_a.print("Ea");
    //E_b.print("Eb");
    //C_a.print("Ca");
    //C_b.print("Cb");

    double E_el = E_electron(N,H,P_a,P_b,F_a,F_b);
    double E_nr = E_nuc_rep(n_atoms,Z,xyz_list);
    double E_tot = E_el + E_nr;
    if (print_energy){
        cout << "Nuclear Repulsion Energy is " << E_nr << " eV" << endl;
        cout << "Electron Energy is " << E_el << " eV" << endl;
        cout << "This molecule in file " << file_name << " has energy " << E_tot << " eV" << endl;
    }

    return E_tot;
}

// derivative of 2 electron integral
mat build_gammaAB_RA(int n_atoms, vector<vector<double>> xyz_list,vector<vec> basis_data,vector<int> atom_list) {
    mat gammaAB_RA(3,n_atoms*n_atoms, fill::zeros);

    for (int coord = 0; coord < 3; coord++) {
        mat gamma_row(n_atoms, n_atoms, fill::zeros);

        for (int A = 0; A < n_atoms; A++) {
            for (int B = 0; B < n_atoms; B++) {

                //Define RA and RB
                vector<double> RA_raw = xyz_list[A];
                vector<double> RB_raw = xyz_list[B];
                vec RA = vec{RA_raw[0], RA_raw[1], RA_raw[2]};
                vec RB = vec{RB_raw[0], RB_raw[1], RB_raw[2]};
                double rA = xyz_list[A][coord];
                double rB = xyz_list[B][coord];

                double g = 0;

                // sum over two-electron integrals over primitive Gaussians
                for (int k = 0; k < 3; k++) {
                    for (int kp = 0; kp < 3; kp++) {
                        for (int l = 0; l < 3; l++) {
                            for (int lp = 0; lp < 3; lp++) {
                                double summand;

                                auto [alpha_k, alpha_kp, alpha_l, alpha_lp, dp_k_sA, dp_kp_sA, dp_l_sB, dp_lp_sB] = gamma_atom_constants(
                                        basis_data, atom_list, RA, RB, A, B, k, kp, l, lp);

                                // Compute 0^0 integral
                                // 6-dimensional 2-center integral over the product of 4 primitive s-type Gaussian functions
                                double sigma_A = 1 / (alpha_k + alpha_kp);
                                double sigma_B = 1 / (alpha_l + alpha_lp);

                                double UA = pow((M_PI * sigma_A), 1.5);
                                double UB = pow((M_PI * sigma_B), 1.5);

                                double V2 = 1 / (sigma_A + sigma_B);
                                double V = sqrt(V2);
                                double diff = norm(RA - RB);
                                double diff2 = pow(diff, 2);
                                double T = V2 * diff2;

                                double oo_RA = 0;
                                if (T == 0){
                                    oo_RA = 0;
                                } else {
                                    oo_RA = UA * UB * (rA - rB)/diff2*(-erf(sqrt(T))/diff+2*V/sqrt(M_PI)*exp(-T));
                                }
                                summand = dp_k_sA * dp_kp_sA * dp_l_sB * dp_lp_sB * oo_RA;
                                g += summand;
                            }
                        }
                    }
                }
                gamma_row(A, B) = 27.2114079527 * g; //convert to eV from atomic units
            }
        }
        gammaAB_RA.row(coord) = gamma_row.as_row();
    }
    return gammaAB_RA;
}

// 1D Analytical Overlap Intergral Derivative
double dI_dxA(double XA, double XB, double alpha, double beta, double lA, double lB){
    return -lA * SxABComp(XA,XB,alpha,beta,lA-1,lB) + 2*alpha * SxABComp(XA,XB,alpha,beta,lA+1,lB);
}

// The derivative of overlap integral
mat build_s_RA(int N, vector<int> basis_atom_nums, vector<tuple<vec,vec,vec,vec,vec>> basis_func_constants){
    // Find all normalization constants for the basis functions
    mat All_Normalization_Constants(3,N,fill::zeros);

    // Run a loop over all your basis functions
    for (int i = 0; i < N; i++){
        // get the normalization constants for the 3 primitives that make up each basis function
        auto[R_center,quantNums,exponents,contraCoeffs,normConsts] = basis_func_constants[i];
        // and save them in an array.
        All_Normalization_Constants.col(i) = normConsts;
    }


    // derivative of Contracted overlap integral S
    mat S_RA(3,N*N, fill::zeros);

    // Express the overlap integral of contracted GTOs into linear
    // combination of primitive GTOs under STO-3G basis set

    for (int coord = 0; coord < 3; coord++) {
        mat S_row(N, N, fill::zeros);

        // Iterate over pairs of basis functions
        for (int mu = 0; mu < N; mu++) {
            for (int nu = 0; nu < N; nu++) {
                // For this given pair of basis functions
                // R_center,quantNums,exponents,contraCoeffs,normConsts for each atom in the pair
                auto [R_mu, lmn_mu, a_mu, d_mu, N_mu] = basis_func_constants[mu];
                auto [R_nu, lmn_nu, a_nu, d_nu, N_nu] = basis_func_constants[nu];

                // # of contracted gaussians
                int K = d_mu.size();
                int L = d_nu.size();

                double S_mu_nu = 0;

                for (int k = 0; k < K; k++) {
                    for (int l = 0; l < L; l++) {
                        double d_k_mu = d_mu(k);
                        double d_l_nu = d_nu(l);
                        double N_k_mu = N_mu(k);
                        double N_l_nu = N_nu(l);
                        double alpha_k_mu = a_mu(k);
                        double alpha_l_nu = a_nu(l);

                        double Xm = R_mu(0);
                        double Ym = R_mu(1);
                        double Zm = R_mu(2);

                        double Xn = R_nu(0);
                        double Yn = R_nu(1);
                        double Zn = R_nu(2);

                        double lm = lmn_mu(0);
                        double mm = lmn_mu(1);
                        double nm = lmn_mu(2);

                        double ln = lmn_nu(0);
                        double mn = lmn_nu(1);
                        double nn = lmn_nu(2);

                        double Ix = SxABComp(Xm, Xn, alpha_k_mu, alpha_l_nu, lm, ln); // Compute Sxkl
                        double Iy = SxABComp(Ym, Yn, alpha_k_mu, alpha_l_nu, mm, mn); // Compute Sykl
                        double Iz = SxABComp(Zm, Zn, alpha_k_mu, alpha_l_nu, nm, nn); // Compute Szkl

                        double S_k_l_RA;

                        if (basis_atom_nums[mu] == basis_atom_nums[nu]){
                            S_k_l_RA = 0;
                        } else {
                            if (coord == 0){ // x
                                double dIdxA = dI_dxA(Xm, Xn, alpha_k_mu, alpha_l_nu, lm, ln);
                                S_k_l_RA = dIdxA * Iy * Iz;
                            } else if (coord == 1){ //y
                                double dIdyA = dI_dxA(Ym, Yn, alpha_k_mu, alpha_l_nu, mm, mn);
                                S_k_l_RA = Ix * dIdyA * Iz;
                            } else if (coord == 2){ //z
                                double dIdzA = dI_dxA(Zm, Zn, alpha_k_mu, alpha_l_nu, nm, nn);
                                S_k_l_RA = Ix * Iy * dIdzA;
                            }
                        }

                        S_mu_nu += d_k_mu * d_l_nu * N_k_mu * N_l_nu * S_k_l_RA;
                    }
                }
                S_row(mu, nu) = S_mu_nu;
            }
        }
        S_RA.row(coord) = S_row.as_row();
    }
    return S_RA;
}

mat build_x(mat p_alpha, mat p_beta, int N, vector<int> basis_atom_nums, vector<int> basis_atom_list, mat IAb){
    mat X(N,N,fill::zeros);
    // inputs basis_atom_nums, basis_orbital_list, IAb, n_atoms, p_tot, Z, gamma, p, s
    // Iterate over pairs of basis functions
    mat P_tot_mu_nu = p_alpha + p_beta;
    for (int mu = 0; mu < N; mu++) {
        for (int nu = 0; nu < N; nu++) {

            double nbetaA = -IAb(2,basis_atom_list[mu]-1);
            double nbetaB = -IAb(2,basis_atom_list[nu]-1);

            double x = (nbetaA+nbetaB)*P_tot_mu_nu(mu,nu);
            X(mu,nu) = x;
        }
    }
    return X;
}

mat build_y(mat P_a, mat P_b, vec p_tot, int N, int n_atoms, vector<int> basis_atom_nums, vec Z){
    mat Y(n_atoms, n_atoms, fill::zeros);
    for (int A = 0; A < n_atoms; A++){
        for (int B = 0; B < n_atoms; B++){
            double sum = 0;

            for (int mu = 0; mu < N; mu++){
                for (int nu = 0; nu < N; nu++){
                    if (basis_atom_nums[mu] == A && basis_atom_nums[nu] == B ){
                        sum += P_a(mu,nu) * P_a(mu,nu) + P_b(mu,nu) * P_b(mu,nu);
                    }
                }
            }
            Y(A,B) = p_tot(A)*p_tot(B) - Z(B)*p_tot(A) - Z(A)*p_tot(B) - sum;
        }
    }
    return Y;


}

// Analytical Gradient/Force of SCF Energy (CNDO/2 Method) - electron term
mat build_E_RA_el(vector<int> basis_atom_nums, mat Xuv, mat Suv_RA, mat YAB, mat gammaAB_RA){
    int n_atoms = YAB.n_rows;
    int N = Xuv.n_rows;
    mat E_RA(3,n_atoms,fill::zeros);
    for (int coord = 0; coord < 3; coord++){
        for (int A = 0; A < n_atoms; A++){
            double term1 = 0;
            double term2 = 0;
            for (int mu = 0; mu < N; mu++){
                for (int nu = 0; nu < N; nu++){
                    if(basis_atom_nums[mu] == A && basis_atom_nums[nu] != A){
                        term1 += Xuv(mu,nu)*Suv_RA(coord,mu*N+nu);
                    }
                }
            }
            for (int B = 0; B < n_atoms; B++){
                if (B != A){
                    term2 += YAB(A,B)*gammaAB_RA(coord,A*n_atoms+B);
                }
            }
            E_RA(coord,A) = term1 + term2;
        }
    }
    return E_RA;
}

// Analytical Gradient/Force of SCF Energy (CNDO/2 Method) - nuclear repulsion term
// derivative of the nuclear repulsion energy
mat build_E_RA_nuc(vec Z, vector<int> basis_atom_nums, vector<vector<double>> xyz_list){
    int n_atoms = Z.n_elem;
    int N = basis_atom_nums.size();
    mat E_RA(3,n_atoms,fill::zeros);
    for (int coord = 0; coord < 3; coord++){
        for (int A = 0; A < n_atoms; A++){
            double sum = 0;
            for (int B = 0; B < n_atoms; B++){
                if (B != A){
                    vector<double> RA_raw = xyz_list[A];
                    vector<double> RB_raw = xyz_list[B];
                    vec RA = vec{RA_raw[0], RA_raw[1], RA_raw[2]};
                    vec RB = vec{RB_raw[0], RB_raw[1], RB_raw[2]};
                    double rA = xyz_list[A][coord];
                    double rB = xyz_list[B][coord];

                    double diff = norm(RA - RB);

                    sum += Z(A)*Z(B)*(rA-rB)/pow(diff,3);
                }
            }
            E_RA(coord,A) = -sum*27.2114079527;
        }
    }
    return E_RA;
}

vector<vector<double>> make_basis_xyz_list(vector<vector<double>> xyz_list,vector<int> basis_atom_nums){
    vector<vector<double>> basis_xyz_list;
    for (int i = 0; i < basis_atom_nums.size(); i++){
        basis_xyz_list.push_back(xyz_list[basis_atom_nums[i]]);
    }
    return basis_xyz_list;
}

vector<int> make_basis_atom_list(vector<int> atom_list,vector<int> basis_atom_nums){
    vector<int> basis_atom_list;
    for (int i = 0; i < basis_atom_nums.size(); i++){
        basis_atom_list.push_back(atom_list[basis_atom_nums[i]]);
    }
    return basis_atom_list;
}

// Analytical Gradient/Force of SCF Energy (CNDO/2 Method) - Full calculation from scratch
mat analytical_E_RA(mat pos, tuple<vector<vec>,vec,vector<int>,vector<int>,string> const_data,
                    bool print_output, bool print_energy){
    //Process inputs
    vector<vec> basis_data = get<0>(const_data);
    vec Z = get<1>(const_data);
    vector<int> atom_list = get<2>(const_data);
    vector<int> basis_atom_nums = get<3>(const_data);
    string file_name = get<4>(const_data);
    int N = basis_atom_nums.size();
    int n_atoms = Z.n_elem;
    vector<vector<double>> xyz_list = Mat2VecVec(n_atoms,pos);
    vector<vector<double>> basis_xyz_list = make_basis_xyz_list(xyz_list,basis_atom_nums);
    vector<int> basis_atom_list = make_basis_atom_list(atom_list,basis_atom_nums);

    // Build Basis Functions
    auto[basis_func_constants,basis_orbital_list] = build_basis(basis_data,N,basis_atom_list,basis_xyz_list);

    // Semi-empirical parameters of ionization energies and electron affinities
    mat IAb = mat("7.176 0 0 0 0 14.051 19.316 25.390 32.272;  0 0 0 0 0 5.572 7.275 9.111 11.080; 9  0 0 0 0 21 25 31 39");

    // Build S, AO basis Contracted overlap integral matrix
    mat S = build_s(N,basis_func_constants);
    mat S_RA = build_s_RA(N,basis_atom_nums,basis_func_constants);

    mat gammaAB = build_g(n_atoms,xyz_list,basis_data,atom_list);
    mat gamma_RA = build_gammaAB_RA(n_atoms,xyz_list,basis_data,atom_list);

    int p = ceil(sum(Z)/2),  q = floor(sum(Z)/2);

    mat H = build_h(N,n_atoms,basis_atom_nums,basis_atom_list,basis_orbital_list,Z,IAb,gammaAB,S);

    auto[P_a,P_b,p_tot,F_a,F_b,C_a,C_b] = CNDO_converge(N, n_atoms, basis_atom_nums,basis_atom_list,S,gammaAB,IAb,basis_orbital_list,Z,p,q);
    double energy  = CNDO_energy(F_a,F_a,N,n_atoms,H,P_a,P_b,Z,xyz_list,file_name, print_energy);

    mat X = build_x(P_a, P_b, N, basis_atom_nums, basis_atom_list, IAb);
    mat Y = build_y(P_a, P_b, p_tot, N, n_atoms, basis_atom_nums, Z);
    mat E_RA_nuc = build_E_RA_nuc(Z,basis_atom_nums,xyz_list);
    mat E_RA_el = build_E_RA_el(basis_atom_nums, X, S_RA, Y, gamma_RA);
    mat E_RA_tot = E_RA_el+E_RA_nuc;

    if (print_output){
        S_RA.print("Suv_RA (A is the center of u)");
        gamma_RA.print("gammaAB_RA");
        E_RA_nuc.print("gradient (Nuclear part)");
        E_RA_el.print("gradient (Electron part)");
        E_RA_tot.print("gradient");
    }

    return E_RA_tot;
}

// From position matrix to energy
double pos2E(mat pos, tuple<vector<vec>,vec,vector<int>,vector<int>,string> const_data,bool print_energy){

    //Process inputs
    vector<vec> basis_data = get<0>(const_data);
    vec Z = get<1>(const_data);
    vector<int> atom_list = get<2>(const_data);
    vector<int> basis_atom_nums = get<3>(const_data);
    string file_name = get<4>(const_data);
    int N = basis_atom_nums.size();
    int n_atoms = Z.n_elem;
    vector<vector<double>> xyz_list = Mat2VecVec(n_atoms,pos);
    vector<vector<double>> basis_xyz_list = make_basis_xyz_list(xyz_list,basis_atom_nums);
    vector<int> basis_atom_list = make_basis_atom_list(atom_list,basis_atom_nums);

    // Build Basis Functions
    auto[basis_func_constants,basis_orbital_list] = build_basis(basis_data,N,basis_atom_list,basis_xyz_list);

    // Semi-empirical parameters of ionization energies and electron affinities
    mat IAb = mat("7.176 0 0 0 0 14.051 19.316 25.390 32.272;  0 0 0 0 0 5.572 7.275 9.111 11.080; 9  0 0 0 0 21 25 31 39");

    // Build S, AO basis Contracted overlap integral matrix
    mat S = build_s(N,basis_func_constants);
    mat S_RA = build_s_RA(N,basis_atom_nums,basis_func_constants);

    mat gammaAB = build_g(n_atoms,xyz_list,basis_data,atom_list);
    mat gamma_RA = build_gammaAB_RA(n_atoms,xyz_list,basis_data,atom_list);

    int p = ceil(sum(Z)/2),  q = floor(sum(Z)/2);

    mat H = build_h(N,n_atoms,basis_atom_nums,basis_atom_list,basis_orbital_list,Z,IAb,gammaAB,S);

    auto[P_a,P_b,p_tot,F_a,F_b,C_a,C_b] = CNDO_converge(N, n_atoms, basis_atom_nums,basis_atom_list,S,gammaAB,IAb,basis_orbital_list,Z,p,q);
    double energy  = CNDO_energy(F_a,F_a,N,n_atoms,H,P_a,P_b,Z,xyz_list,file_name, print_energy);

    return energy;
}

// Golden Section Line Search (from HW 1)
double GoldenSectionSearch(double EnergyConvergenceThreshold,
                           int n,mat F, mat A,
                           double a, double b, double c,tuple<vector<vec>,vec,vector<int>,vector<int>,string> const_data){
    double golden = (3-sqrt(5))/2;
    mat B = A + UnitGrad(n,F)*b;
    mat C = A + UnitGrad(n,F)*c;
    double E_A = pos2E(A,const_data,false);
    double E_B = pos2E(B,const_data,false);
    double E_C = pos2E(C,const_data,false);

    if (abs(E_A-E_B) <= EnergyConvergenceThreshold or abs(E_B-E_C) <= EnergyConvergenceThreshold){
        cout << "Energy converged! *b" << endl;
        return b;
    } else if (abs(E_A-E_C) <= EnergyConvergenceThreshold) {
        cout << "Energy converged! *a" << endl;
        return a;
    } else {
        double x = a + golden*(c-a);
        mat X = A + UnitGrad(n,F)*x;
        double E_X = pos2E(X,const_data,false);
        // if f(x) < f(b)
        if (E_X < E_B){
            // if f(x) < f(b) and a < b < x, new (a,b,c)=(b,x,c)
            if (x > b){
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,b,x,c,const_data);
            } else { // if f(x) < f(b) and a < x < b, new (a,b,c)=(a,x,b)
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,x,b,const_data);
            }
        } else {     // if f(x) > f(b)
            // if f(x) > f(b) and a < b < x, new (a,b,c)=(a,b,x)
            if (x > b){
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,b,x,const_data);
            } else { // if f(x) > f(b) and a < x < b, new (a,b,c)=(x,b,c)
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,x,b,c,const_data);
            }
        }
    }
}

// 1D Steepest Descent (from HW 1)
mat steepest_descent(mat A, tuple<vector<vec>,vec,vector<int>,vector<int>,string> const_data, int& i){

    i += 1; // This is the i-th iteration
    double h = 0.0001;
    double s = 0.3; // initial step
    double ForceConvergenceThreshold = 0.01;
    double EnergyConvergenceThreshold = 1e-8;
    int n = A.n_cols;

    // 1. Input initial point guess (Point A)
    double a = 0;
    double E_A = pos2E(A,const_data,false);

    // 2.1. Compute the gradient, take a step in the opposite direction, check for convergence

    // Compute the gradient
    bool useAnalytical = false;
    mat F;
    if (useAnalytical){
        F = analytical_E_RA(A,const_data,false,false);
    }
    // add else statement for central difference

    // Print New point & grad
    if (!(i==0)){
        cout << "Start golden section search" << endl;
        A.print("new_point");
        cout << "current energy: " << pos2E(A,const_data,false) << endl;
        //mat Cen_F = cF(n,A,ats,h);
        //Cen_F.print("Central Difference Force");
    }


    // check for convergence
    if (GradNorm(n,F) <= ForceConvergenceThreshold){
        cout << "Force converged!" << endl;
        return A;
    } else {
        // take a step in the opposite direction
        double b = s;
        mat B = A + UnitGrad(n,F)*b;
        double E_B = pos2E(B,const_data,false);

        // 2.2. Search for Point B in the opposite direction of the gradient (E(B) < E(A))
        while (!(E_B<E_A)){
            b /= 2;
            B = A + UnitGrad(n,F)*b;
            E_B = pos2E(B,const_data,false);
        }

        // 2.3. Search for Point C in the opposite direction of the gradient (E(C) > E(B))
        double c = 10;
        mat C = A + UnitGrad(n,F)*c;
        double E_C = pos2E(C,const_data,false);
        while (!(E_C > E_B)){
            c *= 2;
            C = A + UnitGrad(n,F)*c;
            E_C = pos2E(C,const_data,false);
        }

        // 2.5. If can find both B and C, use golden section search
        // for the minimum in the opposite direction of the gradient,
        double step_size = GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,b,c,const_data);
        A = A + UnitGrad(n,F)*step_size;

        // return to step 2.1
        return steepest_descent(A,const_data,i);
    }
}

int main(int argc, char* argv[]) {
    // Program inputs
    if (argc !=2){
        printf("Usage: hw3 <filename>, for example hw3 example.txt\n");
        return EXIT_FAILURE;
    }
    string file_name = argv[1];
    //string file_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW3/sample_input/C2H4.txt";

    string H_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW5/basis/H_STO3G.txt";
    string C_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW5/basis/C_STO3G.txt";
    string N_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW5/basis/N_STO3G.txt";
    string O_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW5/basis/O_STO3G.txt";
    string F_path_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW5/basis/F_STO3G.txt";

    // Read input file
    auto[n_atoms,Z,a,b,xyz_list,
         atom_list,basis_xyz_list,
         basis_atom_list,basis_atom_nums] = read_input_file(file_name);

    // Evaluate the number of basis functions, N from the molecular formula, Ca Hb , where
    // the relation is N = 4a + b. Your matrices, such as S, H, X, etc, will be N × N , so this will
    // enable you to define them.
    int N = 4*a+b;

    vector<vec> basis_data = ProcessBasisData(H_path_name,C_path_name,
                                              N_path_name, O_path_name, F_path_name);

    tuple<vector<vec>,vec,vector<int>,
            vector<int>,string> const_data = make_tuple(basis_data, Z, atom_list,
                                                        basis_atom_nums, file_name);

    // Starting coordinates
    mat A = VecVec2Mat(n_atoms,xyz_list);
    mat F = analytical_E_RA(A, const_data, true, true);

    bool optimize = true;
    if (optimize){
        cout << "Start steepest descent with golden section line search using analytical force" << endl;
        int it = 0;
        A = steepest_descent(A,const_data,it);

        cout << "Total iterations: " << it << endl;
        double FinalE = pos2E(A,const_data,false);
        A.print("This is A");
        cout << "Final energy: " << FinalE << endl;
        cout << "Optimized structure: " << endl;
        for (int i = 0; i < n_atoms; ++i) {
            cout << atom_list[i] <<  '(' << A(0,i) << ',' << A(1,i) << ',' << A(2,i) << ')' << endl;
        }
    }
}
