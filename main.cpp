#include <iostream>
#include <iomanip>
#include <armadillo>
#include <random>
#include <cmath>

using namespace std;
using namespace arma;

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
    double J = 0.2;       // Strength of interaction
    double beta = 0.1;    // Temperature
    double s = 0.5;       // Size of spin
    // Initial state
    /*
    vector<double> state = vector<double>(N);  // Initializing the state vector
    for(int i=0; i<N; i++)    state(i) = s; // All spins point up. Could choose something else, of course.
    */

    mat state = mat(N,N);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)    state(i,j) = s; // All spins point up
    }

    // Should I have a vector for the energies? Or use bitwise operations? Probably easier?
    // If I use bitwise operations, I could store one state in one element of a vector. But that is a lot of
    // states. So no.


    double prob;
    double drawn;
    double energy_new;
    double energy_change;
    int P_change;
    double energy_curr = -z*J*s*s*N;


    // Equilibration procedure (Not sure I need this if the initial has all spins pointing up)

    std::default_random_engine generator;                       // I asked the internet, and it replied
    std::uniform_real_distribution<double> distribution(0,1);

    for(int i=0; i<1000; i++)
    {  // Traversing through the spins non-randomly
        cout << "i = " << i << endl;
        for(int j=0; j<N;j++)
        {
            for(int k=0; k<N; k++)
            {   // Know thy neighbours
                // But seriously, there should be a simpler way to do this...
                // Double check at least

                if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                {
                    if(k==0)         energy_change = J*state(0,k)*(state(N-1,k) + state(0,N-1) + state(0,1) + state(1,k));
                    else if(k==N-1)  energy_change = J*state(0,k)*(state(N-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                    else             energy_change = J*state(0,k)*(state(N-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                }
                if(k==0)
                {
                    if(j==N-1)                     energy_change = J*state(j,0)*(state(j-1,0) + state(j,N-1) + state(0,0) + state(0,1));
                    if(j!=0 && j!=N-1)             energy_change = J*state(j,0)*(state(j-1,0) + state(j,N-1) + state(j+1,0) + state(j,1));
                }
                if(j==N-1 && k!=0 && k!=N-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                if(k==N-1 && j!=0 && j!=N-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));

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

    int runs = 10000;
    double msq_av_bin1 = 0;
    double msq_av_bin2 = 0;
    double msq_av_bin3 = 0;
    double m = 0;

    for(int i=0; i<runs; i++)
    {  // Traversing through the spins non-randomly
        for(int j=0; j<N;j++)
        {
            for(int k=0; k<N; k++)
            {   // Know thy neighbours
                // But seriously, there should be a simpler way to do this...
                // Double check at least
                if(j==0)            // Plus sign in front as state(j,k) --> -state(j,k).
                {
                    if(k==0)         energy_change = J*state(0,k)*(state(N-1,k) + state(0,N-1) + state(0,1) + state(1,k));
                    else if(k==N-1)  energy_change = J*state(0,k)*(state(N-1,k) + state(0,k-1) + state(1,k) + state(0,0));
                    else             energy_change = J*state(0,k)*(state(N-1,k) + state(0,k-1) + state(1,k) + state(0,k+1));
                }
                if(k==0)
                {
                    if(j==N-1)                     energy_change = J*state(j,0)*(state(j-1,0) + state(j,N-1) + state(0,0) + state(0,1));
                    if(j!=0 && j!=N-1)             energy_change = J*state(j,0)*(state(j-1,0) + state(j,N-1) + state(j+1,0) + state(j,k+1));
                }
                if(j==N-1 && k!=0 && k!=N-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(0,k) + state(j,k+1));
                if(k==N-1 && j!=0 && j!=N-1)       energy_change = J*state(0,k)*(state(j-1,k) + state(j,k-1) + state(j+1,k) + state(j,0));

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
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)    m += state(i,j);
        }

        msq_av_bin1 += m*m;  // should I maybe take the absolute value instead?
        cout << "m = " << m << "; m**2 = " << msq_av_bin1 << endl;
        m = 0; // Resetting for the next run
    }

    msq_av_bin1 = msq_av_bin1/runs; // Dividing by the number of configurations.

    cout << msq_av_bin1 << endl;    // Crazy high...



} // End of main


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
