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
    double J = 0.2;       // Strength of interaction
    double beta = 0;    // Temperature/ Initializing beta
    int Nbetas = 10;     // Number of different temperatures probed
    double betamin = 0.01;   // Lowest temperature
    double betamax = 100;  // Highest temperature
    double s = 0.5;       // Size of spin
    // Initial state
    /*
    vector<double> state = vector<double>(N);  // Initializing the state vector
    for(int i=0; i<N; i++)    state(i) = s; // All spins point up. Could choose something else, of course.
    */

    // Opening a file to print to
    ofstream printFile;
    string filenamePrefix = "testrun_discard";
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_IsingMC.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    printFile.open(filename);


    printFile << N << " " << s << " " << bins << endl;  // A line of data about the system at the top

    //run_Metropolis(L, dim, z, mcsteps, bins, J, beta, s, printFile);



    vec betas = linspace(betamin, betamax, Nbetas);

    // Where should I put this?
    std::default_random_engine generator;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(0,1);

    // Should I have a vector for the energies? Or use bitwise operations? Probably easier?
    // If I use bitwise operations, I could store one state in one element of a vector. But that is a lot of
    // states. So no.

    double prob;
    double drawn;
    double energy_new;
    double energy_change;
    double energy_curr = -z*J*s*s*N;


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
            for(int j=0; j<L;j++)
            {
                for(int k=0; k<L; k++)
                {   // Know thy neighbours
                    // But seriously, there should be a simpler way to do this...
                    // Double check at least

                    if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                    {
                        if(k==0)         energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,L-1) + state(0,1) + state(1,k));
                        else if(k==L-1)  energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                        else             energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                    }
                    if(k==0)
                    {
                        if(j==L-1)                     energy_change = 2*J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(0,0) + state(0,1));
                        if(j!=0 && j!=L-1)             energy_change = 2*J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(j+1,0) + state(j,1));
                    }
                    if(j==L-1 && k!=0 && k!=L-1)       energy_change = 2*J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                    if(k==L-1 && j!=0 && j!=L-1)       energy_change = 2*J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));
                    if(j!=0 && j!=L-1 && k!=0 && k!=L-1)    energy_change = 2*J*state(j,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,k+1));

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
        }

        //cout << "Done with the equilibration part (now I just have to wait and see how well it actually works...)" << endl;

        // Monte Carlo steps

        //int runs = 10000;
        double msq_av_bin = 0;
        double m = 0;

        for(int l=0; l<bins; l++)
        {
            for(int i=0; i<mcsteps; i++)
            {  // Traversing through the spins non-randomly
                for(int j=0; j<L;j++)
                {
                    for(int k=0; k<L; k++)
                    {   // Know thy neighbours
                        // But seriously, there should be a simpler way to do this...
                        // Double check at least
                        // Lots of steps to avoid segmentation faults
                        // + enforcing boundary conditions
                        if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                        {
                            if(k==0)         energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,L-1) + state(0,1) + state(1,k));
                            else if(k==L-1)  energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                            else             energy_change = 2*J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                        }
                        if(k==0)
                        {
                            if(j==L-1)                     energy_change = 2*J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(0,0) + state(0,1));
                            if(j!=0 && j!=L-1)             energy_change = 2*J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(j+1,0) + state(j,k+1));
                        }
                        if(j==L-1 && k!=0 && k!=L-1)       energy_change = 2*J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                        if(k==L-1 && j!=0 && j!=L-1)       energy_change = 2*J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));
                        if(j!=0 && j!=L-1 && k!=0 && k!=L-1)    energy_change = 2*J*state(j,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,k+1));


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
                    } // End loop over k
                } // End loop over j. All lattice points have been traversed over
                //Commands to find quantity. Magnetization?
                for(int i=0; i<L; i++)
                {
                    for(int j=0; j<L; j++)    m += state(i,j);
                }  // Finding the desired quantity

                msq_av_bin += m*m/(N*N);  // should I maybe take the absolute value instead?
                //cout << "m = " << m << "; m**2 = " << msq_av_bin1 << endl;
                m = 0; // Resetting for the next run
            }  // End over Monte Carlo steps

            msq_av_bin = msq_av_bin/mcsteps; // Dividing by the number of configurations.

            /*
            // Basically, just print commands ... Should write to file instead.
            cout << "Bin no. : " << l << endl;
            cout << "Average magnetization^2: " <<  msq_av_bin << endl;
            cout << "Total magnetization^2: " << msq_av_bin*N*N << endl;
            cout << "RMS Total magnetization: " << sqrt(msq_av_bin*N*N) << endl;
            cout << "Total magnetization^2/s^2: " << msq_av_bin*N*N/(s*s) << endl;
            cout << "RMS Total magnetization*s: " << sqrt(msq_av_bin*N*N*s*s) << endl;
            cout << "Average magnetization^2/s^2: " <<  msq_av_bin/(s*s) << endl;
            cout << "Number of spins: " << N << endl << endl;
            */

            printFile << beta << " " << msq_av_bin << endl;

            msq_av_bin = 0;
        }  // End of loop over bins

    }  // End of loop over betas
    delete filename;
    printFile.close(); // Closing the file

} // End of main

/*
void run_Metropolis(int L, int dim, int z, int mcsteps, int bins, double J, double beta, double s, ofstream printFile)
{
    int N = pow(L, dim);  // Number of sites

    mat state = mat(L,L);
    for(int i=0; i<L; i++)
    {
        for(int j=0; j<L; j++)    state(i,j) = s; // All spins point up
    }

    // Should I have a vector for the energies? Or use bitwise operations? Probably easier?
    // If I use bitwise operations, I could store one state in one element of a vector. But that is a lot of
    // states. So no.


    double prob;
    double drawn;
    double energy_new;
    double energy_change;
    int P_change;         // ?
    double energy_curr = -z*J*s*s*N;


    // Equilibration procedure (Not sure I need this if the initial has all spins pointing up)

    std::default_random_engine generator;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(0,1);

    for(int i=0; i<1000; i++)
    {  // Traversing through the spins non-randomly
        cout << "i = " << i << endl;
        for(int j=0; j<L;j++)
        {
            for(int k=0; k<L; k++)
            {   // Know thy neighbours
                // But seriously, there should be a simpler way to do this...
                // Double check at least

                if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                {
                    if(k==0)         energy_change = J*state(0,k)*(state(L-1,k) + state(0,L-1) + state(0,1) + state(1,k));
                    else if(k==L-1)  energy_change = J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                    else             energy_change = J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                }
                if(k==0)
                {
                    if(j==L-1)                     energy_change = J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(0,0) + state(0,1));
                    if(j!=0 && j!=L-1)             energy_change = J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(j+1,0) + state(j,1));
                }
                if(j==L-1 && k!=0 && k!=L-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                if(k==L-1 && j!=0 && j!=L-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));

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
    }

    cout << "Done with the equilibration part (now I just have to wait and see how well it actually works...)" << endl;

    // Monte Carlo steps

    //int runs = 10000;
    double msq_av_bin = 0;
    double m = 0;

    for(int l=0; l<bins; l++)
    {
        for(int i=0; i<mcsteps; i++)
        {  // Traversing through the spins non-randomly
            for(int j=0; j<L;j++)
            {
                for(int k=0; k<L; k++)
                {   // Know thy neighbours
                    // But seriously, there should be a simpler way to do this...
                    // Double check at least
                    if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                    {
                        if(k==0)         energy_change = J*state(0,k)*(state(L-1,k) + state(0,L-1) + state(0,1) + state(1,k));
                        else if(k==L-1)  energy_change = J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                        else             energy_change = J*state(0,k)*(state(L-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                    }
                    if(k==0)
                    {
                        if(j==L-1)                     energy_change = J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(0,0) + state(0,1));
                        if(j!=0 && j!=L-1)             energy_change = J*state(j,0)*(state(j-1,0) + state(j,L-1) + state(j+1,0) + state(j,k+1));
                    }
                    if(j==L-1 && k!=0 && k!=L-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                    if(k==L-1 && j!=0 && j!=L-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));

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
                } // End loop over k
            } // End loop over j. All lattice points have been traversed over
            //Commands to find quantity. Magnetization?
            for(int i=0; i<L; i++)
            {
                for(int j=0; j<L; j++)    m += state(i,j);
            }

            msq_av_bin += m*m/(N*N);  // should I maybe take the absolute value instead?
            //cout << "m = " << m << "; m**2 = " << msq_av_bin1 << endl;
            m = 0; // Resetting for the next run
        }

        msq_av_bin = msq_av_bin/mcsteps; // Dividing by the number of configurations.

        // Basically, just print commands ... Should write to file instead.
        cout << "Bin no. : " << l << endl;
        cout << "Average magnetization^2: " <<  msq_av_bin << endl;
        cout << "Total magnetization^2: " << msq_av_bin*N*N << endl;
        cout << "RMS Total magnetization: " << sqrt(msq_av_bin*N*N) << endl;
        cout << "Total magnetization^2/s^2: " << msq_av_bin*N*N/(s*s) << endl;
        cout << "RMS Total magnetization*s: " << sqrt(msq_av_bin*N*N*s*s) << endl;
        cout << "Average magnetization^2/s^2: " <<  msq_av_bin/(s*s) << endl;
        cout << "Number of spins: " << N << endl << endl;

        printFile << beta << " " << msq_av_bin << endl;

        msq_av_bin = 0;
    }

} // End of Monte Carlo function
*/



/*
// No, no, no, there is a much simpler way.
double energy(int N, double J, vector<double> state)  // Periodic boundary conditions
{   // But how to determine which ones are neighbours?
    double E = 0;
    for(int i=0; i<(N-1); i++)    E -= J*state(i)*state(i+1);
    E -= J*state(N-1)*state(0);
}

// Have E_change and E_prev and change by

*/
