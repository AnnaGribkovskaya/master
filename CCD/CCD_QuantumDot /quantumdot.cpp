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
    m_ReferenceEnergyHF = HF_Energy;
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
    //cout << eigval[1] << endl;

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

void QuantumDot::SetUpTwoBobyMatrixForHartreeFock(){

    int NumberOfStates = m_shells.size();
    m_twoBodyElementsInHF = create4dArray(NumberOfStates, NumberOfStates, NumberOfStates, NumberOfStates);

    for(int a = 0; a < NumberOfStates; a++){

        for (int b = 0; b < NumberOfStates; b++){

            for(int c = 0; c< NumberOfStates; c++){

                for (int d = 0; d< NumberOfStates; d++){

                    for(int i = 0; i < NumberOfStates; i++) {
                        QuantumState quantum_state_alpha = m_shells.at(i);
                        int alpha_sm = quantum_state_alpha.sm();

                        for(int j = 0; j < NumberOfStates; j++) {
                                QuantumState quantum_state_beta = m_shells.at(j);
                                int beta_sm = quantum_state_beta.sm();

                                for(int k = 0; k < NumberOfStates; k++) {
                                    QuantumState quantum_state_gama = m_shells.at(k);
                                    int gama_sm = quantum_state_gama.sm();

                                    for(int l = 0; l < NumberOfStates; l++) {
                                        QuantumState quantum_state_delta = m_shells.at(l);
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
                                            m_twoBodyElementsInHF[a][b][c][d] += TBME*m_C(a,i)*m_C(b,k)*m_C(c,j)*m_C(d,l);
                                    }
                                }
                        }
                    }

                }

            }
        }

    }
}


void QuantumDot::SetUpInitialAmplitudes(){
    int NumberOfStates = m_shells.size();
    m_twoBodyElementsInHF = create4dArray(NumberOfStates, NumberOfStates, NumberOfStates, NumberOfStates);
    m_InitialAmplitudes = create4dArray(NumberOfStates, NumberOfStates, NumberOfStates, NumberOfStates);
    //applyHartreeFockMethod();
    //SetUpTwoBobyMatrixForHartreeFock();
    for(int i = 0; i < NumberOfParticles; i++) {
        //QuantumState quantum_state_alpha = m_shells.at(i);
        //int alpha_sm = quantum_state_alpha.sm();

         for(int j = 0; j < NumberOfParticles; j++) {
             //QuantumState quantum_state_beta = m_shells.at(j);
             //int beta_sm = quantum_state_beta.sm();

             for(int a = NumberOfParticles; a < NumberOfStates; a++) {
                 //QuantumState quantum_state_gama = m_shells.at(k);
                 //int gama_sm = quantum_state_gama.sm();

                 for(int b = NumberOfParticles; b < NumberOfStates; b++) {
                     //QuantumState quantum_state_delta = m_shells.at(l);
                     //int delta_sm = quantum_state_delta.sm();
                     m_InitialAmplitudes[a][b][i][j] = (m_twoBodyElementsInHF[a][b][i][j] - m_twoBodyElementsInHF[a][b][j][i])/(eigval[i]+eigval[j]-eigval[a]-eigval[b]);



                  }
              }
           }
      }

}


void QuantumDot:: ComputeCorrelationEnergy(){
    m_CorrelationEnergy = 0.0;
    int NumberOfStates = m_shells.size();
    for(int i = 0; i < NumberOfParticles; i++) {
         for(int j = 0; j < NumberOfParticles; j++) {
             for(int a = NumberOfParticles; a < NumberOfStates; a++) {
                 for(int b = NumberOfParticles; b < NumberOfStates; b++) {

                     m_CorrelationEnergy += 0.25*m_Amplitudes[a][b][i][j]*m_twoBodyElementsInHF[i][j][a][b];




                  }
              }
           }
      }

}

void QuantumDot::ComputeEpsilon(){
    int NumberOfStates = m_shells.size();
    for(int i = 0; i < NumberOfParticles; i++) {
         for(int j = 0; j < NumberOfParticles; j++) {
             for(int a = NumberOfParticles; a < NumberOfStates; a++) {
                 for(int b = NumberOfParticles; b < NumberOfStates; b++) {
                     m_Epsilon[a][b][i][j] = 1/(eigval[i]+eigval[j]-eigval[a]-eigval[b]);
                  }
              }
           }
      }


}

void QuantumDot:: ComputeAmplitudes(){
    int NumberOfStates = m_shells.size();
    ComputeEpsilon();
    m_Amplitudes = create4dArray(NumberOfStates, NumberOfStates, NumberOfStates, NumberOfStates);
    m_Amplitudes_previous = create4dArray(NumberOfStates, NumberOfStates, NumberOfStates, NumberOfStates);
    //cout << m_Epsilon<< endl;
    for(int i = 0; i < NumberOfParticles; i++) {
         for(int j = 0; j < NumberOfParticles; j++) {
             for(int a = NumberOfParticles; a < NumberOfStates; a++) {
                 for(int b = NumberOfParticles; b < NumberOfStates; b++) {
                     double first_term = 0.0;
                     double second_term = 0.0;
                     double third_term = 0.0;
                     double forth_term = 0.0;
                     double last_term1 = 0.0;
                     double last_term2 = 0.0;
                     double last_term3 = 0.0;
                     double last_term4 = 0.0;

                     first_term = m_twoBodyElementsInHF[a][b][i][j];

                             for(int c = NumberOfParticles; c < NumberOfStates; c++) {
                                 for(int d = NumberOfParticles; d < NumberOfStates; d++) {
                                     second_term+=m_twoBodyElementsInHF[a][b][c][d]*m_Amplitudes_previous[c][d][i][j];
                                 }
                             }

                             for(int k = 0; k < NumberOfParticles; k++) {
                                 for(int l = 0; l < NumberOfParticles; l++) {
                                     third_term+=m_twoBodyElementsInHF[k][l][i][j]*m_Amplitudes_previous[a][b][k][l];
                                 }
                             }

                             for(int k = 0; k < NumberOfParticles; k++) {
                                 for(int c = NumberOfParticles; c < NumberOfStates; c++) {
                                     forth_term+=m_twoBodyElementsInHF[k][b][c][j]*m_Amplitudes_previous[a][c][i][k]-m_twoBodyElementsInHF[k][a][c][j]*m_Amplitudes_previous[b][c][i][k]-
                                             m_twoBodyElementsInHF[k][b][c][i]*m_Amplitudes_previous[a][c][j][k]-m_twoBodyElementsInHF[k][a][c][i]*m_Amplitudes_previous[b][c][j][k];
                                 }
                             }

                             for(int k = 0; k < NumberOfParticles; k++) {
                                 for(int l = 0; l < NumberOfParticles; l++) {
                                     for(int c = NumberOfParticles; c < NumberOfStates; c++) {
                                         for(int d = NumberOfParticles; d < NumberOfStates; d++) {
                                             last_term1+=m_twoBodyElementsInHF[k][l][c][d]*m_Amplitudes_previous[c][d][i][j]*m_Amplitudes_previous[a][b][k][l];
                                             last_term2+=m_twoBodyElementsInHF[k][l][c][d]*(m_Amplitudes_previous[a][c][i][k]*m_Amplitudes_previous[b][d][j][l]-m_Amplitudes_previous[a][c][j][k]*m_Amplitudes_previous[b][d][i][l]);
                                             last_term3+=m_twoBodyElementsInHF[k][l][c][d]*(m_Amplitudes_previous[d][c][i][k]*m_Amplitudes_previous[a][b][l][j]-m_Amplitudes_previous[d][c][j][k]*m_Amplitudes_previous[a][b][l][i]);
                                             last_term4+=m_twoBodyElementsInHF[k][l][c][d]*(m_Amplitudes_previous[a][c][l][k]*m_Amplitudes_previous[d][b][i][j]-m_Amplitudes_previous[b][c][l][k]*m_Amplitudes_previous[d][a][i][j]);
                                         }
                                     }
                                 }
                             }
                     m_Amplitudes[a][b][i][j]=m_Epsilon[a][b][i][j]*(first_term+0.5*second_term+0.5*third_term+forth_term+0.25*last_term1+last_term2-0.5*last_term3-0.5*last_term4);
                     cout << m_Amplitudes[a][b][i][j] << endl;
                  }
              }
           }
      }

}


void QuantumDot:: applyCoupledClusterDoubles(){
    applyHartreeFockMethod();
    SetUpTwoBobyMatrixForHartreeFock();
    SetUpInitialAmplitudes();
    m_Amplitudes=m_InitialAmplitudes;
    //cout<< m_Amplitudes << endl;
    ComputeCorrelationEnergy();
    cout << "===================================================================" << endl;
    cout << "Reference Energy:  " <<  m_ReferenceEnergyHF << endl;
    cout << "Correlation Energy  " << m_CorrelationEnergy << endl;
    m_Amplitudes_previous = m_InitialAmplitudes;
    ComputeAmplitudes();
    cout<< m_Amplitudes<< endl;
    ComputeCorrelationEnergy();
    cout << "Correlation Energy Next  " << m_CorrelationEnergy << endl;
    double energy_difference = 1000;
    double energy_next = 0;
    double energy_prev = 0;
    double epsilon = 10e-8;
    int i = 0;
    while (epsilon < energy_difference && i < 1000){
        ComputeCorrelationEnergy();
        energy_prev =  m_CorrelationEnergy;
        m_Amplitudes = m_Amplitudes_previous;
        ComputeAmplitudes();
        ComputeCorrelationEnergy();
        energy_next = m_CorrelationEnergy;
        energy_difference = abs(energy_next - energy_prev);
        i++;
    }
    ComputeCorrelationEnergy();
    cout << "===================================================================" << endl;
    cout << "Reference Energy:  " <<  m_ReferenceEnergyHF << endl;
    cout << "Correlation Energy  " << m_CorrelationEnergy << endl;
    cout << "Ground State Energy  " << m_ReferenceEnergyHF + m_CorrelationEnergy << endl;
    cout << "===================================================================" << endl;
    cout<< "Last computed energy"<< m_ReferenceEnergyHF + energy_next << endl;

}
