#include <iostream>
#include <iomanip>
#include <armadillo>
#include <random>
#include <cmath>
#include <string>

using namespace std;
using namespace arma;
using std::ofstream; using std::string;

void run_Metropolis(int L, int dim, int z, int mcsteps, int bins, double J, double beta, double s, ofstream printFile);
double energy(int N, double J, vector<double> state);
double energychanged(double energy);


int main()
{
    // Initialization //
    // System properties
    int L = 15;           // Size of system
    int dim = 2;          // Dimension of system
    int z = 4;            // Number of neighbours
    int N = pow(L, dim);  // Number of sites   
    int mcsteps = 1000;   // Number of Monte Carlo steps
    int bins = 100;       // Number of bins
    double J = 1;         // Strength of interaction
    double beta;          // Temperature/ Initializing beta
    int Nbetas = 100;     // Number of different temperatures probed
    double betamin = 0;   // Lowest temperature
    double betamax = 5;  // Highest temperature
    double s = 0.5;       // Size of spin


    // Opening a file to print to
    ofstream printFile;
    string filenamePrefix = "beta0to5_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_cverrorattempt";
    //string filenamePrefix = "test";
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_IsingMC.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    printFile.open(filename);

    ofstream allFile;
    char *filename2 = new char[1000];
    sprintf(filename2, "%s_IsingMC_all.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    allFile.open(filename2);


    printFile << N << " " << s << " " << bins << endl;  // A line of data about the system at the top
    allFile << N << " " << s << " " << bins << " " << Nbetas << " " << mcsteps << " " << betamin << " " << betamax << endl;  // A line of data about the system at the top


    vec betas = linspace(betamin, betamax, Nbetas);             // For simulating different temperatures

    // Magnetization
    vec msq_av_bin = zeros(bins);
    mat binholder_m = mat(bins,mcsteps);                        // Holds the values for each MC step.

    // Energy
    vec av_E = zeros(bins);
    mat binholder_E = mat(bins, mcsteps);

    // Energy squared
    vec av_Esq = zeros(bins);
    mat binholder_Esq = mat(bins, mcsteps);

    // Heat capacity
    vec av_c = zeros(bins);



    // Where should I put this?
    std::default_random_engine generator;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(0,1);

    std::default_random_engine generator2;
    std::uniform_int_distribution<int> distribution2(0,L-1);

    // Quantities for determining system state
    int x;
    int y;
    int xp1;
    int xm1;
    int yp1;
    int ym1;
    double prob;
    double drawn;
    double energy_new;
    double energy_change;


    // Relevant quantities
    double m_average;    // Well, <m^2>, m is magnetization per site.
    double E_average;    // The energy for the total system.
    double Esq_average;  // Squared energy for the total system.
    double cv;           // Found by E_average and Esq_average.

    for(int n=0; n<Nbetas; n++)
    {
        // Picking the beta
        beta = betas(n);

        // Setting up the intial state
        mat state = mat(L,L);
        for(int i=0; i<L; i++)
        {
            for(int j=0; j<L; j++)    state(i,j) = s; // All spins point up
        }
        double energy_curr = -(z/2)*J*s*s*N;

        // Equilibration procedure


        for(int i=0; i<1000; i++)
        {  // Traversing through the spins randomly
            for(int j=0; j<N; j++)
            {
                x = distribution2(generator2);
                y = distribution2(generator2);
                // Periodic boundary conditions
                if(x==(L-1))                    xp1 = 0;
                else                            xp1 = x+1;
                if(x==0)                        xm1 = L-1;
                else                            xm1 = x-1;
                if(y==(L-1))                    yp1 = 0;
                else                            yp1 = y+1;
                if(y==0)                        ym1 = L-1;
                else                            ym1 = y-1;

                energy_change = 2*J*state(x,y)*(state(xm1,y) + state(x,ym1) + state(xp1,y) + state(x,yp1));
                energy_new = energy_curr + energy_change;

                if(energy_new <= energy_curr)
                {
                    state(x,y) = -state(x,y); // Flip spin if the result is a state with lower energy.
                    energy_curr = energy_new; // Update energy
                }
                else
                {
                    prob = exp(-beta*(energy_new-energy_curr));
                    drawn = distribution(generator);
                    if(drawn<prob)
                    {

                        state(x,y) = -state(x,y);  // Flip spin.
                        energy_curr = energy_new;  // Update energy
                    }
                }
            }
        }

        //cout << "Done with the equilibration part (now I just have to wait and see how well it actually works...)" << endl;

        // Monte Carlo steps

        //int runs = 10000;
        double m = 0;

        //double counter = 0;
        for(int l=0; l<bins; l++)
        {   // Loop over the bins
            msq_av_bin(l) = 0;
            av_E(l)       = 0;
            av_Esq(l)     = 0;
            av_c(l)       = 0;
            for(int i=0; i<mcsteps; i++)
            {   // Loop over the Monte Carlo steps
                // Traversing through the spins randomly
                for(int j=0; j<N;j++)
                {
                    x = distribution2(generator2);
                    y = distribution2(generator2);
                    // Periodic boundary conditions
                    if(x==(L-1))                    xp1 = 0;
                    else                            xp1 = x+1;
                    if(x==0)                        xm1 = L-1;
                    else                            xm1 = x-1;
                    if(y==(L-1))                    yp1 = 0;
                    else                            yp1 = y+1;
                    if(y==0)                        ym1 = L-1;
                    else                            ym1 = y-1;

                    energy_change = 2*J*state(x,y)*(state(xm1,y) + state(x,ym1) + state(xp1,y) + state(x,yp1));
                    energy_new = energy_curr + energy_change;

                    if(energy_new <= energy_curr)
                    { // If the new energy is lower than the energy of the old, we commit the change
                        state(x,y) = -state(x,y); // Flip spin if the result is a state with lower energy.
                        energy_curr = energy_new; // Update energy
                    }
                    else
                    { // If the new energy is larger than the energy of the old, we leave it to chance
                        prob = exp(-beta*(energy_new-energy_curr));
                        drawn = distribution(generator);
                        if(drawn<prob)
                        {
                            state(x,y) = -state(x,y);  // Flip spin.
                            energy_curr = energy_new;  // Update energy
                        }
                    } // End if-tests of energy
                } // End loop over j. N lattice points have been traversed over
                //Commands to find quantity
                m = 0;    // Resetting for this run
                for(int a=0; a<L; a++)
                {
                    for(int b=0; b<L; b++)    m += state(a,b);
                }  // Finding the desired quantity

                double yo = m*m/(N*N); // Man gave name to all the quantities. In the beginning. In the beginning.
                msq_av_bin(l) += yo;   // Accumulating the average for bin l
                binholder_m(l,i)= yo;  // Save each value of m grouped by bin.

                av_E(l) += energy_curr;
                binholder_E(l,i) = energy_curr;

                av_Esq(l) += energy_curr*energy_curr;
                binholder_Esq(l,i) = energy_curr*energy_curr;

                allFile << m << " " << energy_curr << endl;

            }  // End of Monte Carlo steps. Index i out.
            // Dividing by the number of configurations to get the average for bin l
            msq_av_bin(l) = msq_av_bin(l)/mcsteps;
            av_E(l)       = av_E(l)/mcsteps;
            av_Esq(l)     = av_Esq(l)/mcsteps;
            av_c(l)       = beta*beta*(av_Esq(l)-av_E(l)*av_E(l));

            //cout << av_E(l) << endl;

        }  // End of loop over bins. Index l out.

        // Magnetization
        // Calculating the average magnetization over all bins.
        m_average = 0;
        for(int i=0; i<bins; i++)    m_average += msq_av_bin(i);
        m_average = m_average/bins;

        // Use this:
        // What the others seem to be suggesting (Formula from Kosovan and Sega)
        // This is quite similar to the variance...
        double blockvariance = 0;
        // Accumulating
        for(int i=0; i<bins; i++)    blockvariance += (msq_av_bin(i)-m_average)*(msq_av_bin(i)-m_average);
        // Dividing
        blockvariance = blockvariance/(bins*(bins-1));
        double block_stdv = sqrt(blockvariance);

        // Energy and heat capacity
        // Energy //
        E_average = 0;
        for(int i=0; i<bins; i++)    E_average += av_E(i);
        E_average = E_average/bins;

        // Error in the energy //
        double E_stdv = 0;
        for(int i=0; i<bins; i++)    E_stdv += (av_E(i)-E_average)*(av_E(i)-E_average);
        E_stdv = sqrt(E_stdv/(bins*(bins-1)));

        // Energy squared //
        Esq_average = 0;
        for(int i=0; i<bins; i++)    Esq_average += av_Esq(i);
        Esq_average = Esq_average/bins;

        // Error in the energy //
        double Esq_stdv = 0;
        for(int i=0; i<bins; i++)    Esq_stdv += (av_Esq(i)-Esq_average)*(av_Esq(i)-Esq_average);
        Esq_stdv = sqrt(Esq_stdv/(bins*(bins-1)));

        cv = beta*beta*(Esq_average-E_average*E_average);

        // Is this right?
        double cv_average  = 0;
        for(int i=0; i<bins; i++)    cv_average += av_c(i);
        cv_average = cv_average/bins;

        // Should I use cv_average or cv in the error estimation?
        double cv_stdv = 0;
        for(int i=0; i<bins; i++)    cv_stdv += (av_c(i)-cv_average)*(av_c(i)-cv_average);
        cv_stdv = sqrt(cv_stdv/(bins*(bins-1)));

        // Printing to file
        printFile << beta << " " << m_average << " " << block_stdv << " " << cv << " " << E_average << " " << Esq_average << " " << E_stdv << " " << Esq_stdv << " " << cv_stdv << " " << cv_average << endl;


    }  // End of loop over betas. Index n out.
    delete filename;
    printFile.close(); // Closing the file

    delete filename2;
    allFile.close();

} // End of main

