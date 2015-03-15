/*===========================================================================
 * dynamic_nets_direct.cc
 *---------------------------------------------------------------------------
 * Simulate dynamical processes over a network structure that varies over 
 * time in accordance with direct and indirect (delayed) co-occurance data
 * in space.
 * 
 * Last Update: 15th May 2012
 * Version:     2.0
 *===========================================================================*/

#include <dynamic_nets_direct.h>
#include <netevo.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace lemon;
using namespace netevo;

/**
 * Global of the nodes in the network (each is an ant).
 */
vector<Node> nodes;

/**
 * Print program usage.
 */
void printUsage (void) {
   cout << "SI Spread Over Dynamic Networks Driven By Data (Version 1.0)" << endl;
   cout << "Usage: dynNetDirect FILENAME SIZE SI_PROB DECAY_RATE ANT RUNS LEN PREFIX" << endl;
   cout << "  FILENAME:   Interaction data file." << endl;
   cout << "  SIZE:       Number of ants in data file." << endl;
   cout << "  SI_PROB:    S->I transition probability." << endl;
   cout << "  ANT:        Ant to start infected (-1 = run for all)." << endl;
   cout << "  RUNS:       Number of randomised runs." << endl;
   cout << "  LEN:        Timesteps per simulation." << endl;
   cout << "  TIMESTEP:   Length of a time step." << endl;
   cout << "  OUT_FREQ:   Output frequency (timesteps)." << endl;
   cout << "  PREFIX:     Prefix for output files." << endl;
}

/** 
 * Calculates the weight that edge should have given a delayed crossing.
 * t is the time period, a is the rate of decay. Use of an exponential
 * function ensures that result is in range (0, 1) for t >= 0
 */
double calcWeight (double t, double a) {
   return exp(-a*t);
}

/** 
 * SI Dynamics.
 * Uses the dynamic network from data to influence spread. We use the 
 * calcWeight() function to calculate a decay of the S->I probability
 * given a particular delay since last crossing.
 */
class SIMap : public NodeDynamic {
protected:
   double m_probSI;
   double m_decayRate;
   DynamicNet &m_net;
   double m_ts;
public:   
   SIMap (double probSI, double decayRate, DynamicNet &net, double ts) : m_probSI(probSI), 
      m_decayRate(decayRate), m_net(net), m_ts(ts) { }
   string getName () { return "SIMap"; }
   int getStates () { return 1; } // (0 = Suseptible, 1 = Infected)
   void setDefaultParams (Node v, System &sys) { }
   
   void fn (Node v, System &sys, const State &x, State &dx, const double t) {
      int i;
      double tt;
      int vID = sys.stateID(v);
      double prob, rndNum, crossing;
      
      tt = m_ts * t;
      
      // Only consider uninfected nodes
      if (x[vID] == 0.0) {
         // Search through all possible neighbours to see if infected
         for (i=0; i<m_net.getSize(); ++i) {
            if (i != vID && x[i] == 1.0) {
            	crossing = m_net.checkInteraction(i, vID, tt, tt+m_ts);
               if (crossing != -1.0) {
                  if (sys.rnd() <= m_probSI) {
                     // An infection has occured, stop searching any further
                     dx[vID] = 1.0;
                     // Update the infected time
                     m_net.setInfectedTime(vID, tt);
                     return;
                  }
               }
            }
         }
      }
      
      // Nothing has changed.
      dx[vID] = x[vID];
   }
};

/**
 * Run simulations for a particular ant and output to a given prefix.
 * This will output to file the results for a given ant. Each run is in
 * a separate file.
 */
void doRuns (System &sys, DynamicNet &dynNet, int ant, int runs, double simLen, double ts, int outFreq, const char *prefix) {
   int i, j, k;
   char buf[1000];
   
   // Generate a random initial state for the simulation
   State initial = State(sys.totalStates(), 0.0);
   initial[ant] = 1.0;
   dynNet.setInfectedTime(ant, 0.0);

   // Create a simulator for mapping dynamics
   SimulateMap simMap;
   
   // Output the data to a vector for later writing to file
   vector<double> tOut;
   vector<State>  xOut;
   SimObserverToVectors vectorObserver(xOut, tOut);

   // We don't need to log changes so use the default change logger that does nothing
   ChangeLog nullLogger;

   // Create and open the file for output
   sprintf(buf, "%sANT-%i.txt", prefix, ant+1);
   ofstream outFile;
   outFile.open(buf);
   
   // Simulate for the required number of times.
   for (i=0; i<runs; ++i) {
      
      // Copy the initial state and clear output vectors
      State initialCopy = initial;
      tOut.clear();
      xOut.clear();
      
      // Simulate the dynamics for our initial state and using the observer and logger
      simMap.simulate(sys, simLen, initialCopy, vectorObserver, nullLogger);
      
      // Save the simulation results to file.
      for (j=0; j<tOut.size(); ++j) {
         if (j%outFreq == 0 || j == (tOut.size() - 1)) {
            State curState = xOut[j];
            outFile << (i+1) << "," << tOut[j] * ts << "," << (int)(curState[0]);
            for (k=1; k<curState.size(); ++k) {
               outFile << "," << (int)(curState[k]);
            }
            outFile << endl;
         }
      }
   }
   
   // Finish writing to file
   outFile.close();
}

/** 
 * Main function.
 */
int main (int argc, const char **argv) {
   double probSI, decayRate, simLen, ts;
   int num, ant, runs, i, outFreq;
   const char *netFile, *prefix;
   
   // Check that there is a correct number of arguments.
   if (argc < 9 || argc > 11) {
      printUsage();
      return 1;
   }
   
   // Gather all the command line arguments (convert if necessary).
   netFile = argv[1];
   num = atoi(argv[2]);
   probSI = atof(argv[3]);
   //decayRate = atof(argv[4]);
   ant = atoi(argv[4]);
   runs = atoi(argv[5]);
   simLen = atof(argv[6]);
   ts = atof(argv[7]);
   outFreq = atoi(argv[8]);
   if (argc == 10) { prefix = argv[9]; }
   
   // Create a dynamic network structure used by the dynamics.
   // Must provide size of network and file name.
   DynamicNet net = DynamicNet(num, netFile);
   
   // Create the system.
   System sys;
   
   // Load the required dynamics.
   SIMap vDyn(probSI, decayRate, net, ts);
   sys.addNodeDynamic(&vDyn);
   
   // Add the nodes with the SI dynamics.
   for (i=0; i<num; ++i) {
      nodes.push_back(sys.addNode("SIMap"));
   }
   
   // Run the simulations for an individual ant or all ants.
   if (ant == -1) {
      for (i=0; i<num; ++i) {
         doRuns(sys, net, i, runs, simLen, ts, outFreq, prefix);
      }
   }
   else if (ant > 0 && ant <= num) {
      doRuns(sys, net, ant-1, runs, simLen, ts, outFreq, prefix);
   }
   else {
      cerr << "Error: incorrect ant number specified." << endl;
      return 1;
   }
   
   // Ended successfully.
   return 0;
}
