#include "quantumdot.h"
#include <iostream>


using namespace std;

bool myComparison(const pair<int,double> &a,const pair<int,double> &b)
{
       return a.second<b.second;
}

QuantumDot::QuantumDot(int EnergyCutOff, double h_omega, int ParticlesNumber){
    setUpStatesPolarSorted(EnergyCutOff, h_omega, ParticlesNumber);
}


void QuantumDot::setUpStatesCartesian(int EnergyCutOff) {
    QuantumState m_q_state;
    int nx, ny;
    for (int i = 1; i < EnergyCutOff + 1; i++){
        for (int j = 0; j < i; j++){
            nx = j;
            ny = i - nx -1;
            m_q_state.set(nx, ny, m_sm, m_s);
            m_shells.push_back(m_q_state);
            m_q_state.flipSpin();
            m_shells.push_back(m_q_state);
         }
    }
}

void QuantumDot::setUpStatesPolarSorted(int EnergyCutOff, double h_omega, int ParticlesNumber) {
    homega = h_omega;
    NumberOfParticles = ParticlesNumber;
    m_EnergyCutOff = EnergyCutOff;
    QuantumState m_q_state;
    vector<int> oddShells;
    int n, m;

    //loop to find all odd shell numbers including energyCutOff shell
    for (int i = 1; i <= EnergyCutOff; i++){
        if (i % 2 != 0){
            oddShells.push_back(i);
        }
    }
    for (int index = 0; index < oddShells.size(); ++index){
       //positive m
       for (int j=0; j<=EnergyCutOff - oddShells[index]; j++){
           n = index;
           m = j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
       }
       //negative m
       for (int j=1; j<=EnergyCutOff - oddShells[index]; j++){
           n = index;
           m = -1*j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
       }
    }

    //now sort m_shells vector
    pair<int,double> mapping;
    vector<pair<int,double>> vector_to_sort;
    for(int i = 0; i < m_shells.size(); i++) {
        QuantumState quantum_state = m_shells.at(i);
        double EnergyOfState = homega*((double)2.0*quantum_state.n() + abs(quantum_state.m()) + 1.0);
        mapping = make_pair(i, EnergyOfState);
        vector_to_sort.push_back(mapping);
    }

    sort(vector_to_sort.begin(),vector_to_sort.end(),myComparison);
    vector<QuantumState> sorted_states;
    for(int i = 0; i < vector_to_sort.size(); i++) {
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
        m_shells.at(vector_to_sort.at(i).first).flipSpin();
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
    }
    m_shells = sorted_states;
    int dim = m_shells.size();
    CalculateNonIntEnergy();
    m_twoBodyElements = create4dArray(dim, dim, dim, dim);
    m_CCD_t_old = create4dArray(dim, dim, dim, dim);
    m_CCD_t = create4dArray(dim, dim, dim, dim);
    fillTwoBodyElements();

}

double**** QuantumDot::create4dArray(int dim1, int dim2, int dim3, int dim4) {

    double**** h = new double***[dim1];

    for(int i=0; i < dim1; i++) {
        h[i] = new double**[dim2];
        for (int j=0; j < dim2; j++) {
            h[i][j] = new double*[dim3];
            for(int k = 0; k < dim3; k++) {
                h[i][j][k] = new double[dim4];
            }
        }
    }
    for(int i=0; i < dim1; i++) {
        for (int j=0; j < dim2; j++) {
            for(int k = 0; k < dim3; k++) {
                for(int l = 0; l < dim4; l++) {
                    h[i][j][k][l] = 0;
                }

            }
        }
    }
    return h;
}

void QuantumDot::fillTwoBodyElements(){
    int NumberOfStates = m_shells.size();

    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();

        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();

            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();

                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();

                    if ((alpha_sm == beta_sm && gama_sm == delta_sm)){
                        m_twoBodyElements[i][k][j][l] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, beta_n, beta_m,  delta_n, delta_m);
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm)){
                        m_twoBodyElements[i][k][l][j] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, delta_n, delta_m, beta_n, beta_m);
                    }
                }
            }
        }
    }
}

void QuantumDot::setCoefficientMatrix(arma::mat CoefficientMatrix){
     m_C = CoefficientMatrix;
}

arma::mat QuantumDot::computeDensityMatrix(){
    int NumberOfStates = m_shells.size();
    arma::mat DensityMatrix(NumberOfStates, NumberOfStates);

    for(int k = 0; k < NumberOfStates; k++) {
        for(int l = 0; l < NumberOfStates; l++) {
            double sum = 0.0;
            for (int i=0; i < NumberOfParticles; i++) {
                    sum += m_C(k,i)*m_C(l,i);
                    DensityMatrix(k, l) = sum;
            }
        }
    }
    return DensityMatrix;
}

void QuantumDot::CalculateNonIntEnergy(){
    int NumberOfStates = m_shells.size();
    m_HOEnergies.zeros(NumberOfStates,NumberOfStates);
    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state = m_shells.at(i);
        m_HOEnergies(i, i) = (2.0*(double)quantum_state.n() + (double)abs(quantum_state.m()) + 1.0)*homega;
    }
}

void QuantumDot::computeHFmatrix(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    m_HF.zeros(NumberOfStates,NumberOfStates);
    double FockElement = 0;

    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();

        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();

            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();

                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();
                    double TBME = 0.0;
                    double tbme1 = 0.0;
                    double tbme2 = 0.0;

                    if ((alpha_sm == beta_sm && gama_sm == delta_sm)){
                        tbme1 = m_twoBodyElements[i][k][j][l];
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm)){
                        tbme2 = m_twoBodyElements[i][k][l][j];
                    }
                    TBME = tbme1 - tbme2;
                    FockElement += DensityMatrix(k,l)*TBME;
                }
            }
            if (i == j) {
                m_HF(i, i) += m_HOEnergies(i, i);
            }
            m_HF(i, j) += FockElement;
            FockElement = 0.0;
        }
    }
}

double QuantumDot::computeHartreeFockEnergyDifference(){
    return ((arma::accu(abs(eigval - eigval_previous)))/(double)m_shells.size());
}

void QuantumDot::computeHartreeFockEnergy(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;
    double selfConsistentFIeldIterations = 0.0;
    double ExchangePart = 0.0;
    double SingleParticleEnergies = 0.0;

    for(int f = 0; f < FermiLevel; f++){
        SingleParticleEnergies += eigval(f);
    }

    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();

        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();

            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();

                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();

                    double TBME = 0.0;
                    double tbme1 = 0.0;
                    double tbme2 = 0.0;
                    if ((alpha_sm == beta_sm) && (gama_sm == delta_sm)){
                       tbme1 = m_twoBodyElements[i][k][j][l];
                    }
                    if ((alpha_sm == delta_sm) && (gama_sm == beta_sm)){
                       tbme2 = m_twoBodyElements[i][k][l][j];
                    }
                    TBME = tbme1 - tbme2;
                    selfConsistentFIeldIterations = DensityMatrix(i,j)*DensityMatrix(k,l)*TBME;
                    ExchangePart += selfConsistentFIeldIterations;
                }
            }
        }
    }
    double HF_Energy = SingleParticleEnergies - 0.5*ExchangePart;
    // Uncoment for debug
    //cout << "SPEnergies " << SingleParticleEnergies << endl;
    //cout << "Exchange " << ExchangePart << endl;
    cout << "===================================================================" << endl;
    cout << setprecision(12);
    cout << "Num of electrons = " << NumberOfParticles << endl;
    cout << "Num of shells = " << m_EnergyCutOff << endl;
    cout << "Omega = " << homega << endl;
    cout << "Total energy " << HF_Energy << endl;
    writeToFile(HF_Energy, NumberOfParticles, m_EnergyCutOff, homega);
}

void QuantumDot::applyHartreeFockMethod(){
    int NumberOfStates = m_shells.size();
    arma::mat C(NumberOfStates, NumberOfStates);

    C.eye();
    setCoefficientMatrix(C);
    double difference = 10; //dummy value to handle first iteration
    double epsilon = 10e-8;

    eigval_previous.zeros(NumberOfStates);
    int i = 0;
    while (epsilon < difference && i < 1000){
        arma::mat x_DensityMatrix = computeDensityMatrix();
        computeHFmatrix(x_DensityMatrix);
        arma::eig_sym(eigval, eigvec, m_HF);
        setCoefficientMatrix(eigvec);
        difference = computeHartreeFockEnergyDifference();
        eigval_previous = eigval;
        i++;

    }

    arma::mat y_DensityMatrix = computeDensityMatrix();
    m_numOfIterations = i;
    computeHartreeFockEnergy(y_DensityMatrix);
    cout << "Number of iterations " << i << endl;
}

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << "n = " << quantum_state.n() << endl;
        cout << "m = " <<quantum_state.m() << endl;
        cout << "sm = " <<quantum_state.sm() << endl;
        cout << "s = " <<quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

void QuantumDot::getQuantumDotStatesNumber(){
    cout << "Number of available states of system is " << m_shells.size() << endl;
    CalculateNonIntEnergy();
    fillTwoBodyElements();
    setUpFmatrix();
    computeInitialCCDCorrEnergy();
    applyCCDMethod();

}

void QuantumDot::writeToFile(double HF_Energy, int NumberOfParticles, int m_EnergyCutOff, double homega){
    ofstream ofile;
    ofile.open(ResultsFile, ios::app);
    ofile << setprecision(12);
    ofile << "===============================" << endl;
    ofile << "Num of electrons = " << NumberOfParticles << endl;
    ofile << "Num of shells = " << m_EnergyCutOff << endl;
    ofile << "Omega = " << homega << endl;
    ofile << "Total energy " << HF_Energy << endl;
    ofile << "Iterations " << m_numOfIterations << endl;
    ofile << "Eigenvalues: " << endl;
    for (int i = 0; i < eigval.size(); i++){
        ofile << eigval(i) << "     "<< m_HOEnergies(i, i) << endl;
    }
    ofile.close();
}

double QuantumDot::computeInitialCCDCorrEnergy(){
    double CorrelationEnergy = 0.0;
    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;
    for(int i = 0; i < FermiLevel; i++){
        for(int j = 0; j < FermiLevel; j++){
            for(int a = FermiLevel; a < NumberOfStates; a++){
                for(int b = FermiLevel; b < NumberOfStates; b++){
                    CorrelationEnergy +=  (m_twoBodyElements[a][b][i][j]-m_twoBodyElements[a][b][j][i])*
                            (m_twoBodyElements[i][j][a][b]-m_twoBodyElements[i][j][b][a])/( m_HOEnergies(i, i) +  m_HOEnergies(j, j) -  m_HOEnergies(a, a) -  m_HOEnergies(b, b));

                }
            }
        }
    }
    cout << "CORR NEW EN   " << 0.25*CorrelationEnergy << endl;
    return 0.25*CorrelationEnergy;
}

void QuantumDot::computeInitialCCDAmplitudes(){
    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;
    for(int i = 0; i < FermiLevel; i++){
        for(int j = 0; j < FermiLevel; j++){
            for(int a = FermiLevel; a < NumberOfStates; a++){
                for(int b = FermiLevel; b < NumberOfStates; b++){

                    m_CCD_t_old[a][b][i][j] = (m_twoBodyElements[a][b][i][j]-m_twoBodyElements[a][b][j][i])/(m_HOEnergies(i, i) +  m_HOEnergies(j, j) -  m_HOEnergies(a, a) -  m_HOEnergies(b, b));

                }
            }
        }
    }
}

double QuantumDot::computeCCDCorrEnergy(){
    double CorrelationEnergy = 0.0;
    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;
    for(int i = 0; i < FermiLevel; i++){
        for(int j = 0; j < FermiLevel; j++){
            for(int a = FermiLevel; a < NumberOfStates; a++){
                for(int b = FermiLevel; b < NumberOfStates; b++){
                    CorrelationEnergy += (m_twoBodyElements[i][j][a][b]-m_twoBodyElements[i][j][b][a])*m_CCD_t[a][b][i][j];
                }
            }
        }
    }
    return 0.25*CorrelationEnergy;
}

void QuantumDot::computeCCDAmplitudes(){
    double EnergyDenominator;
    int FermiLevel = NumberOfParticles;
    int NumberOfStates = m_shells.size();
    for(int j = 0; j < FermiLevel; j++){
        for(int i = 0; i < FermiLevel; i++){
            for(int b = FermiLevel; b < NumberOfStates; b++){
                for(int a = FermiLevel; a < NumberOfStates; a++){
                    EnergyDenominator = m_HOEnergies(i, i) +  m_HOEnergies(j, j) -  m_HOEnergies(a, a) -  m_HOEnergies(b, b);

                    m_CCD_t[a][b][i][j] = (m_twoBodyElements[a][b][i][j]-m_twoBodyElements[a][b][j][i]);

/*
                    for(int c = FermiLevel; c < NumberOfStates; c++){
                        if( c!=b){
                            m_CCD_t[a][b][i][j] += (m_twoBodyElements[a][c][i][j]-m_twoBodyElements[a][c][j][i])*m_fm(b,c);
                        }
                        if( c!=a){
                            m_CCD_t[a][b][i][j] += (m_twoBodyElements[c][b][i][j]-m_twoBodyElements[c][b][j][i])*m_fm(a,c);
                        }

                    }
                    for(int k = 0; k < FermiLevel; k++){
                        if( k!=j){
                            m_CCD_t[a][b][i][j] -= (m_twoBodyElements[a][b][i][k]-m_twoBodyElements[a][b][k][i])*m_fm(k,j);
                        }
                        if( k!=i){
                            m_CCD_t[a][b][i][j] -= (m_twoBodyElements[a][b][k][j]-m_twoBodyElements[a][b][j][k])*m_fm(k,i);
                        }

                    }

*/
                    for(int c = FermiLevel; c < NumberOfStates; c++){
                        for(int d = FermiLevel; d < NumberOfStates; d++){
                            m_CCD_t[a][b][i][j] += 0.5*(m_twoBodyElements[a][b][c][d]-m_twoBodyElements[a][b][d][c])*m_CCD_t_old[c][d][i][j];
                        }
                    }

                    for(int k = 0; k < FermiLevel; k++){
                        for(int l = 0; l < FermiLevel; l++){
                            m_CCD_t[a][b][i][j] += 0.5*(m_twoBodyElements[k][l][i][j]-m_twoBodyElements[k][l][j][i])*m_CCD_t_old[a][b][k][l];
                        }
                    }

                    for(int k = 0; k < FermiLevel; k++){
                        for(int l = 0; l < FermiLevel; l++){
                            for(int c = FermiLevel; c < NumberOfStates; c++){
                                for(int d = FermiLevel; d < NumberOfStates; d++){
                                    m_CCD_t[a][b][i][j] += 0.25*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*m_CCD_t_old[c][d][i][j]*m_CCD_t_old[a][b][k][l];
                                    m_CCD_t[a][b][i][j] += (m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[a][c][i][k]*m_CCD_t_old[b][d][j][l] - m_CCD_t_old[a][c][j][k]*m_CCD_t_old[b][d][i][l]);
                                    m_CCD_t[a][b][i][j] -= 0.5*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[d][c][i][k]*m_CCD_t_old[a][b][l][j] - m_CCD_t_old[d][c][j][k]*m_CCD_t_old[a][b][l][i]);
                                    m_CCD_t[a][b][i][j] -= 0.5*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[a][c][l][k]*m_CCD_t_old[d][b][i][j] - m_CCD_t_old[b][c][l][k]*m_CCD_t_old[d][a][i][j]);
                                }
                            }
                        }
                    }


                    for(int k = 0; k < FermiLevel; k++){
                        for(int c = FermiLevel; c < NumberOfStates; c++){
                            m_CCD_t[a][b][i][j] += (m_twoBodyElements[k][b][c][j]-m_twoBodyElements[k][b][j][c])*m_CCD_t_old[a][c][i][k] - (m_twoBodyElements[k][b][c][i]-m_twoBodyElements[k][b][i][c])*m_CCD_t_old[a][c][j][k]  - (m_twoBodyElements[k][a][c][j]-m_twoBodyElements[k][a][j][c])*m_CCD_t_old[b][c][i][k] + (m_twoBodyElements[k][a][c][i]-m_twoBodyElements[k][a][i][c])*m_CCD_t_old[b][c][j][k];
                        }
                    }

                    //EnergyDenominator = m_Epsilon(i) + m_Epsilon(j) - m_Epsilon(a) - m_Epsilon(b);



                    //if (m_CCD_t[a][b][i][j] != 0){
                    //    cout << m_CCD_t[a][b][i][j] << endl;
                    //}

                    m_CCD_t[a][b][i][j] /= EnergyDenominator;
                }
            }
        }
    }
}

void QuantumDot::updateOldAmplitudes(){
    int NumberOfStates = m_shells.size();
    for(int i = 0; i < NumberOfStates; i++){
        for(int j = 0; j < NumberOfStates; j++){
            for(int k = 0; k < NumberOfStates; k++){
                for(int l = 0; l < NumberOfStates; l++){
                    m_CCD_t_old[i][j][k][l] = m_CCD_t[i][j][k][l];
                }
            }
        }
    }
}

void QuantumDot::applyCCDMethod(){
    int NumberOfStates = m_shells.size();
    double EnergyDiff = 100.0;
    double CorrelationEnergy;
    double CorrelationEnergyOld = computeInitialCCDCorrEnergy();
    cout << setprecision(16) << "initial correlation energy(MBPT2) " << CorrelationEnergyOld << endl;
    computeInitialCCDAmplitudes();

    double m_toleranceCCD = 0.000001;
    int FermiLevel = NumberOfParticles;
    int i= 0;
    while (abs(EnergyDiff) > m_toleranceCCD && i<10){
        i++;
        cout << "loop start" << endl;
        //cout << "Ediff   :     " << abs(EnergyDiff) << endl;
        computeCCDAmplitudes();
        cout << "-----------computed  amplitudes -----------------" << endl;

        for(int i = 0; i < FermiLevel; i++){
            for(int j = 0; j < FermiLevel; j++){
                for(int a = FermiLevel; a < NumberOfStates; a++){
                    for(int b = FermiLevel; b < NumberOfStates; b++){


                        if (m_CCD_t[a][b][i][j] != 0.0) {
                        cout << "a b i j " << a << " "<< b<<" " <<i <<" " << j<< " " << m_CCD_t[a][b][i][j] << endl;}


                    }
                }
            }
        }
        CorrelationEnergy = computeCCDCorrEnergy();
        cout <<setprecision(16) << "correlation energy in loop " << CorrelationEnergy << endl;
        EnergyDiff = CorrelationEnergy - CorrelationEnergyOld;
        CorrelationEnergyOld = CorrelationEnergy;
        updateOldAmplitudes();

    }
    cout << setprecision (16) << "final correlation energy CCD " << CorrelationEnergy  << endl;
}

void QuantumDot::computeInitialT1Amplitudes(){

    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;

    m_CCD_t1_old.zeros(NumberOfStates,FermiLevel);
    /*
    for(int i = 0; i < FermiLevel; i++){
        for(int a = FermiLevel; a < NumberOfStates; a++){
                    //m_CCD_t1_old(a,i) = (m_twoBodyElements[a][b][i][j]-m_twoBodyElements[a][b][j][i])/(m_HOEnergies(i, i) +  m_HOEnergies(j, j) -  m_HOEnergies(a, a) -  m_HOEnergies(b, b));
        }
    }
    */
}

void QuantumDot::setUpFmatrix(){
     int FermiLevel = NumberOfParticles;
     int NumberOfStates = m_shells.size();
     m_fm.zeros(NumberOfStates,NumberOfStates);
     cout << "Begin set up f-matrix" << endl;
     for(int j = 0; j < NumberOfStates; j++){
         for(int i = 0; i < NumberOfStates; i++){
             m_fm(i,j) = m_HOEnergies(i,j);
             for(int k = 0; k < FermiLevel; k++){
             m_fm(i,j) +=  m_twoBodyElements[i][k][j][k]-m_twoBodyElements[i][k][k][j];
             }
         }
     }
     cout << "End set up f-matrix" << endl;
    //m_fm.print();

}


void QuantumDot::computeT1Amplitudes(){
    double EnergyDenominator;
    int FermiLevel = NumberOfParticles;
    int NumberOfStates = m_shells.size();
    for(int i = 0; i < FermiLevel; i++){
        for(int a = FermiLevel; a < NumberOfStates; a++){

            EnergyDenominator = m_HOEnergies(i, i) -  m_HOEnergies(a, a);
            m_CCD_t1(a,i) = m_fm(a,i);

            for(int k = 0; k < FermiLevel; k++){
                for(int c = FermiLevel; c < NumberOfStates; c++){

                    m_CCD_t1(a,i) += (m_twoBodyElements[k][a][c][i] - m_twoBodyElements[k][a][i][c])*m_CCD_t1_old(c,k); //L1
                    m_CCD_t1(a,i) += m_fm(k,c)*m_CCD_t_old[a][c][i][k];//L2
                    m_CCD_t1(a,i) -= m_fm(k,c)*m_CCD_t1_old(c,i)*m_CCD_t1_old(a,k);//Q1

                    for(int d = FermiLevel; d < NumberOfStates; d++){

                        m_CCD_t1(a,i) += 0.5*(m_twoBodyElements[k][a][c][d] - m_twoBodyElements[k][a][d][c])*m_CCD_t_old[c][d][k][i];//L3
                        m_CCD_t1(a,i) += (m_twoBodyElements[k][a][c][d] - m_twoBodyElements[k][a][d][c])*m_CCD_t1_old(c,k)*m_CCD_t1_old(d,i);//Q3
                    }
                    for(int l = 0; l < FermiLevel; l++){

                        m_CCD_t1(a,i) -= 0.5*(m_twoBodyElements[k][l][c][i] - m_twoBodyElements[k][l][i][c])*m_CCD_t_old[c][a][k][l];//L4
                        m_CCD_t1(a,i) -= (m_twoBodyElements[k][l][c][i] - m_twoBodyElements[k][l][i][c])*m_CCD_t1_old(c,k)*m_CCD_t1_old(a,l);//Q2
                        for(int d = FermiLevel; d < NumberOfStates; d++){

                            m_CCD_t1(a,i) += (m_twoBodyElements[k][l][c][d] - m_twoBodyElements[k][l][d][c])*m_CCD_t1_old(c,k)*(m_CCD_t_old[d][a][l][i] - m_CCD_t1_old(d,i)*m_CCD_t1_old(a,l));//T12
                            m_CCD_t1(a,i) -= 0.5*(m_twoBodyElements[k][l][c][d] - m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[c][d][k][i]*m_CCD_t1_old(a,l) + m_CCD_t_old[c][a][k][l]*m_CCD_t1_old(d,i));//T34
                        }
                    }
                }
            }
            m_CCD_t1(a,i) /= EnergyDenominator;
        }
    }
}

void QuantumDot::computeT2Amplitudes(){
    double EnergyDenominator;
    int FermiLevel = NumberOfParticles;
    int NumberOfStates = m_shells.size();
    for(int j = 0; j < FermiLevel; j++){
        for(int i = 0; i < FermiLevel; i++){
            for(int b = FermiLevel; b < NumberOfStates; b++){
                for(int a = FermiLevel; a < NumberOfStates; a++){
                    EnergyDenominator = m_HOEnergies(i, i) +  m_HOEnergies(j, j) -  m_HOEnergies(a, a) -  m_HOEnergies(b, b);

                    m_CCD_t[a][b][i][j] = (m_twoBodyElements[a][b][i][j]-m_twoBodyElements[a][b][j][i]);


                    for(int c = FermiLevel; c < NumberOfStates; c++){
                        for(int d = FermiLevel; d < NumberOfStates; d++){
                            m_CCD_t[a][b][i][j] += 0.5*(m_twoBodyElements[a][b][c][d]-m_twoBodyElements[a][b][d][c])*m_CCD_t_old[c][d][i][j];
                        }
                    }

                    for(int k = 0; k < FermiLevel; k++){
                        for(int l = 0; l < FermiLevel; l++){
                            m_CCD_t[a][b][i][j] += 0.5*(m_twoBodyElements[k][l][i][j]-m_twoBodyElements[k][l][j][i])*m_CCD_t_old[a][b][k][l];
                        }
                    }

                    for(int k = 0; k < FermiLevel; k++){
                        for(int l = 0; l < FermiLevel; l++){
                            for(int c = FermiLevel; c < NumberOfStates; c++){
                                for(int d = FermiLevel; d < NumberOfStates; d++){
                                    m_CCD_t[a][b][i][j] += 0.25*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*m_CCD_t_old[c][d][i][j]*m_CCD_t_old[a][b][k][l];
                                    m_CCD_t[a][b][i][j] += (m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[a][c][i][k]*m_CCD_t_old[b][d][j][l] - m_CCD_t_old[a][c][j][k]*m_CCD_t_old[b][d][i][l]);
                                    m_CCD_t[a][b][i][j] -= 0.5*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[d][c][i][k]*m_CCD_t_old[a][b][l][j] - m_CCD_t_old[d][c][j][k]*m_CCD_t_old[a][b][l][i]);
                                    m_CCD_t[a][b][i][j] -= 0.5*(m_twoBodyElements[k][l][c][d]-m_twoBodyElements[k][l][d][c])*(m_CCD_t_old[a][c][l][k]*m_CCD_t_old[d][b][i][j] - m_CCD_t_old[b][c][l][k]*m_CCD_t_old[d][a][i][j]);
                                }
                            }
                        }
                    }


                    for(int k = 0; k < FermiLevel; k++){
                        for(int c = FermiLevel; c < NumberOfStates; c++){
                            m_CCD_t[a][b][i][j] += (m_twoBodyElements[k][b][c][j]-m_twoBodyElements[k][b][j][c])*m_CCD_t_old[a][c][i][k] - (m_twoBodyElements[k][b][c][i]-m_twoBodyElements[k][b][i][c])*m_CCD_t_old[a][c][j][k]  - (m_twoBodyElements[k][a][c][j]-m_twoBodyElements[k][a][j][c])*m_CCD_t_old[b][c][i][k] + (m_twoBodyElements[k][a][c][i]-m_twoBodyElements[k][a][i][c])*m_CCD_t_old[b][c][j][k];
                        }
                    }

                    //Difference for CCSD

                    for(int c = FermiLevel; c < NumberOfStates; c++){
                        m_CCD_t[a][b][i][j] += (m_twoBodyElements[a][b][c][j]-m_twoBodyElements[a][b][j][c])*m_CCD_t1_old(c,i) - (m_twoBodyElements[a][b][c][i]-m_twoBodyElements[a][b][i][c])*m_CCD_t1_old(c,j);
                    }


                    for(int k = 0; k < FermiLevel; k++){
                        for(int c = FermiLevel; c < NumberOfStates; c++){
                            m_CCD_t[a][b][i][j] -= (m_twoBodyElements[k][b][i][j]-m_twoBodyElements[j][b][j][i])*m_CCD_t1_old(a,k) - (m_twoBodyElements[k][a][i][j]-m_twoBodyElements[j][a][j][i])*m_CCD_t1_old(b,k);
                        }
                    }







                    m_CCD_t[a][b][i][j] /= EnergyDenominator;
                }
            }
        }
    }
}
