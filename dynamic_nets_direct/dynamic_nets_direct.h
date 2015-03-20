/*
 * dynamic_nets_direct.h
 *
 * Functions to load and allow for arbitary times to be
 * easily calculated.
 */

#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

/**
 * Dynamic network that uses real data to drive edge weights.
 */
class DynamicNet {
private:
   /** Number of nodes in network. */
   int m_size;
   /**
    * Lists of the nodes and their delayed crossing times. 
    * List of list of pairs; 1st list has element for each node;
    * 2nd list is for each individual crossing; pair holds
    * the <crossing time, time other node was there>.
    */
   vector< pair<double,double> > **states;
   
   vector< double > *infectedTime;
   
   /** Gets the state crossing vector for a given edge. */
   vector< pair<double,double> > & getState (int to, int from) {
      return *states[(m_size * from) + to];
   }
     
public:
   /**
    * Constructor for a dynamic data driven network.
    * Must supply the number of nodes (size) and a filename to extract
    * the crossing data from.
    */
   DynamicNet (int size, string filename) { 
      int i, j, from, to;
      m_size = size;
      states = (vector< pair<double,double> > **)malloc(sizeof(vector< pair<double,double> > *) * size * size);
      for (i = 0; i < (size * size); ++i) {
            states[i] = new vector< pair<double,double> >();
      }
      
      infectedTime = new vector< double >(size, -1.0);

      ifstream infile(filename.c_str());

      vector <string> record;
      string none("NA");

      while (infile) {
         string s;
         if (!getline(infile, s)) break;

         istringstream ss(s);
         record.clear();

         while (ss) {
            string s;
            if (!getline(ss, s, '\t')) break;
            record.push_back(s);
         }

         from = atoi(record[0].c_str()) - 1;
         to = atoi(record[1].c_str()) - 1;
         addUpdate(from, to, atof(record[2].c_str()), atof(record[3].c_str()));
      }

      if (!infile.eof()) {
       cerr << "Could not load file.\n";
      }
   };
   
   /**
    * Destructor for the dynamic data driven network.
    * Frees allocated memory.
    */
   ~DynamicNet () { 
      int i, j;
      for (i = 0; i < m_size; ++i) {
         free(states[i]);
      }
      free(states);
      free(infectedTime);
   };
   
   /**
    * Adds a crossing.
    */
   void addUpdate (int from, int to, double fromTime, double toTime) {
      // Assumes that input data is sorted in accending time. If not then
      getState(from, to).push_back(make_pair(fromTime, toTime));
      getState(to, from).push_back(make_pair(fromTime, toTime));
   };
   
   /**
    * Calculate if an interaction is taking place in a given time interval.
    */
   double checkInteraction (int from, int to, double t_start, double t_end) {
      vector< pair<double,double> >::iterator itr;
      double firstTime;
      double l = 0.0, r = 0.0;
      
      // Check to ensure there has been a crossing
      if (getState(from, to).size() == 0) {
         return -1.0;
      }

      // Find the last crossing before the given time and calculate difference
      for (itr = getState(from, to).begin(); itr < getState(from, to).end(); ++itr) {
         if ( ((*itr).first >= t_start && (*itr).first < t_end) ||
              ((*itr).second >= t_start && (*itr).second < t_end) ) { 
            return 1.0;
         }
      }

      // Crossing is not happening at this time point, so ignore.
      return -1.0;
   };
   
   /** Getter and setter for the infected time. */
   double getInfectedTime (int node) { return infectedTime->at(node); }
   void setInfectedTime (int node, double time) { infectedTime->at(node) = time; }
   
   /** Return the number of nodes in the network. */
   int getSize () { return m_size; }
};

