#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

// This program works by reading in pH data from a plain text file (note that pH values should be added only into one column with n rows).
// and reading that pH data into a vector, here. Equilibrium calculations are performed for each element in the vector (each pH measurement)...
// (see ammonia_concentration function) and the resulting ammonia concentration is written to a separate file. 

// Note the constants below for a given buffer system.

const double initial_pH = 6.45; // subject to change

const double buffer_concentration = 1; // molarity
const double buffer_ratio = 0.63095734448017; // basic:acidic
const double monobasic_concentration = 0.61313682015316;
const double dibasic_concentration = (buffer_concentration - monobasic_concentration);
const double monobasic_pKa = initial_pH - log10(buffer_ratio); // subject to change
double Ka = pow(10,-(monobasic_pKa)); // monobasic Ka
double Kb = 1.778 * pow(10,-5); // ammonia Kb
double Kw = pow(10,-14); // water auto-ionization
const double ammonia_monobasic_K = (Ka * Kb)/Kw; // equilibrium constant for nh3 + h2po4 = nh4(+) + HPO4(2-)


// read in pH data: excel to txt to vector
void read_pHdata(vector<double> &pHdata, istream &input) {
    double pH;
    while(input >> pH){
        pHdata.push_back(pH);
    } // end loop
} // end function

// calculate ammonia concentration for a given pH value
void ammonia_concentration(const vector<double> &pHdata, ostream &output){
    double dx;
    double new_buffer_ratio;
    double ammonia_concentration;
    double numer;
    double denom;

    for (size_t i = 0; i < pHdata.size(); ++i) {
        new_buffer_ratio = pow(10,pHdata.at(i)-monobasic_pKa);
        dx = ((new_buffer_ratio*monobasic_concentration) - dibasic_concentration)/(1+new_buffer_ratio);
        numer = (dx*(dx+dibasic_concentration+(monobasic_concentration*ammonia_monobasic_K)-(dx*ammonia_monobasic_K)));
        denom = (ammonia_monobasic_K*(monobasic_concentration-dx));
        ammonia_concentration = numer/denom;

        output << ammonia_concentration << endl; 
    } // end loop
} // end function

int main() {

    vector<double> pHdata;
    ifstream pHInputStream("pHdata_lab.txt");

    if (!(pHInputStream.is_open())){
        cout << "File doesn't exist." << endl;
        return -1;
    }
    
    read_pHdata(pHdata, pHInputStream); // read in pH data
    pHInputStream.close(); // close stream when done

    ofstream ammoniaOutputStream("ammoniaConcentration_data.txt");
    if (!(ammoniaOutputStream.is_open())) {
        cout << "Output stream doesn't exist" << endl;
        return -1;
    }

    ammonia_concentration(pHdata, ammoniaOutputStream);
    ammoniaOutputStream.close();
    
} // end function
