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
    double betamax = 0.1;  // Highest temperature
    double s = 0.5;       // Size of spin
    // Initial state
    /*
    vector<double> state = vector<double>(N);  // Initializing the state vector
    for(int i=0; i<N; i++)    state(i) = s; // All spins point up. Could choose something else, of course.
    */

    // Opening a file to print to
    ofstream printFile;
    string filenamePrefix = "beta0to0p1_Nbetas100_J1_mcsteps1000_bins100";
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_IsingMC.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    printFile.open(filename);

    ofstream allFile;
    char *filename2 = new char[1000];
    sprintf(filename2, "%s_IsingMC_all.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    allFile.open(filename2);


    printFile << N << " " << s << " " << bins << endl;  // A line of data about the system at the top

    //run_Metropolis(L, dim, z, mcsteps, bins, J, beta, s, printFile);

    vec betas = linspace(betamin, betamax, Nbetas);             // For simulate different temperatures
    mat binholder_m = mat(bins,mcsteps);                           // Holds the values for each MC step.
                                                                //To be deleted after each bin

    // Where should I put this?
    std::default_random_engine generator;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(0,1);

    std::default_random_engine generator2;
    std::uniform_real_distribution<int> distribution2(0,L-1);

    // Should I have a vector for the energies? Or use bitwise operations? Probably easier?
    // If I use bitwise operations, I could store one state in one element of a vector. But that is a lot of
    // states. So no.

    // Quantities for determining system state
    double x;
    double y;
    double xp1;
    double xm1;
    double yp1;
    double ym1;
    double prob;
    double drawn;
    double energy_new;
    double energy_change;
    double energy_curr = -z*J*s*s*N;

    //
    double m_average;


    // For error analysis/ correlation functions
    double covariance_term = 0;
    double variance = 0;
    double sigma_msq = 0;
    double autocorrelation_time = 0;

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

        // Equilibration procedure (Not sure I need this if the initial has all spins pointing up)


        for(int i=0; i<1000; i++)
        {  // Traversing through the spins non-randomly
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
                    state(j,k) = -state(j,k); // Flip spin if the result is a state with lower energy.
                    energy_curr = energy_new; // Update energy
                }
                else
                {
                    prob = exp(-beta*(energy_new-energy_curr));
                    drawn = distribution(generator);
                    if(drawn<prob)
                    {

                        state(j,k) = -state(j,k);  // Flip spin.
                        energy_curr = energy_new;  // Update energy
                    }
                }
            }
        }

        //cout << "Done with the equilibration part (now I just have to wait and see how well it actually works...)" << endl;

        // Monte Carlo steps

        //int runs = 10000;
        vec msq_av_bin = zeros(bins);
        double m = 0;

        //double counter = 0;
        for(int l=0; l<bins; l++)
        {   // Loop over the bins
            msq_av_bin(l) = 0;

            for(int i=0; i<mcsteps; i++)
            {   // Loop over the Monte Carlo steps
                // Traversing through the spins non-randomly
                for(int j=0; j<L;j++)
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

                    energy_new = energy_curr + energy_change;
                    if(energy_new <= energy_curr)
                    { // If the new energy is lower than the energy of the old, we commit the change
                        state(j,k) = -state(j,k); // Flip spin if the result is a state with lower energy.
                        energy_curr = energy_new; // Update energy
                    }
                    else
                    { // If the new energy is larger than the energy of the old, we leave it to chance
                        prob = exp(-beta*(energy_new-energy_curr));
                        drawn = distribution(generator);
                        if(drawn<prob)
                        {
                            state(j,k) = -state(j,k);  // Flip spin.
                            energy_curr = energy_new;  // Update energy
                        }
                    } // End if-tests of energy
                } // End loop over j. All lattice points have been traversed over
                //Commands to find quantity
                m = 0;    // Resetting for this run
                for(int a=0; a<L; a++)
                {
                    for(int b=0; b<L; b++)    m += state(a,b);
                }  // Finding the desired quantity

                double yo = m*m/(N*N); // Man gave name to all the quantities. In the beginning. In the beginning.
                msq_av_bin(l) += yo;   // Accumulating the average for bin l
                binholder_m(l,i)= yo;  // Save each value of m grouped by bin.
            }  // End of Monte Carlo steps. Index i out.
            // Dividing by the number of configurations to get the average for bin l
            msq_av_bin(l) = msq_av_bin(l)/mcsteps;

        }  // End of loop over bins. Index l out.

        // Calculating the average magnetization over all bins.
        m_average = 0;
        for(int i=0; i<bins; i++)    m_average += msq_av_bin(i);
        m_average = m_average/bins;

        // Variance, according to MHJ
        variance = 0;
        // Loop for accumulating terms for the variance
        for(int i=0; i<bins; i++)
        {
            for(int k=0; k<mcsteps; k++)    variance += (binholder_m(i,k)-m_average)*(binholder_m(i,k)-m_average);
        }
        // Dividing to get the variance
        variance = variance/(bins*mcsteps);


        // Finding the covariance term in sigma_m, according to MHJ:
        covariance_term = 0;
        // Loop for accumulating terms for the covariance term
        for(int i=0; i<bins; i++)
        {
            for(int k=0; k<mcsteps; k++)
            {
                for(int p=0; p<k; p++)
                {
                    //double hi = (binholder_m(i,k)-m_average)*(binholder_m(i,p)-m_average);
                    covariance_term += (binholder_m(i,k)-m_average)*(binholder_m(i,p)-m_average);
                    /*if(hi<0)
                    {
                        // cout << "Negative contribution to  the covariance term" << endl;
                        counter++;
                        //cout << "This is the " << counter << "th time" << endl;
                    }*/
                }
            }
        }
        // Multiplying to get the covariance term
        covariance_term = covariance_term*2/(bins*mcsteps*mcsteps);

        // Sigma_m^2, according to MHJ
        sigma_msq = variance/mcsteps + covariance_term;

        // Autocorrelation time, according to MHJ
        autocorrelation_time = 0;
        for(int i=0; i<bins; i++)
        {
            for(int d=0; d<(bins-1); d++) // Kind of a weird limit
            {
                for(int k=0; k<(bins-d); k++)    autocorrelation_time += (binholder_m(i,k)-m_average)*(binholder_m(k+d)-m_average);
            }
        }
        autocorrelation_time = autocorrelation_time*2/(bins*mcsteps*variance);
        autocorrelation_time += 1;
        // Hardcoding to avoid infinities in correlation-free experiment
        if(sigma_msq==0 && variance==0 && covariance_term==0)    autocorrelation_time = 1; // Correlation free experiment: autocorrelation_time=1
        // Get infinity when variance=0, but when covariance=0, we have a correlation-free experiment.

        // What the others seem to be suggesting (Formula from Kosovan and Sega)
        // This is quite similar to the variance...
        double blockvariance = 0;
        // Accumulating
        for(int i=0; i<bins; i++)    blockvariance += (msq_av_bin(i)-m_average)*(msq_av_bin(i)-m_average);
        // Dividing
        blockvariance = blockvariance/(bins*(bins-1));


        printFile << beta << " " << m_average << " " << sigma_msq << " " << variance << " " << covariance_term << " " << autocorrelation_time << " " << blockvariance  << endl;

        //cout << "Number of negative contributions to the covariance term:" << counter << endl;
        //counter = 0;

    }  // End of loop over betas. Index n out.
    delete filename;
    printFile.close(); // Closing the file

    delete filename2;
    allFile.close();

} // End of main

