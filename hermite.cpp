//#include "polynomial.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


#include <map>
#include <vector>

using namespace std;

typedef map<double, vector<double> >    Table;
Table table;

bool hermite_read(char * filename)
{
    ifstream file(filename);
    if (!file.good())
    {
        cerr << "Cannot open file " << filename << "!\n";
        return false;
    }
    double x, y;
    vector<double> ys;
    int i = 1;
    for (string line; getline(file, line); ++i)
    {
        // cout << "line " << i;
        if (line == "") continue;
        istringstream iss(line);
        if (!(iss >> x))
        {
            cerr << "Error: Cannot read x value in line " << i << "!\n";
            table.clear();
            return false;
        }
        int j = 1;
        while (iss >> y)
        {
            ys.push_back(y);
            ++j;
        }
        if (!iss.eof())
        {
            cerr << "Error: Cannot read y value " << j << " in line " << i << "!\n";
            table.clear();
            return false;
        }
        table[x] = move(ys);
        // cout << " ok" << endl;
    }
    
    return true;
}

void hermite_show()
{
    for (auto & kvp : table)
    {
        
        double         const & x  = kvp.first;
        vector<double> const & ys = kvp.second;
        cout << setw(3) << x << " : ";
        for (auto const & y : ys)
        {
            cout << setw(3) << y;
        }
        cout << '\n';
    }
    cout << endl;
}

void hermite_ipol()
{
    vector<vector<double>> yss;
    vector<double> ys;
    for (auto const & kvp : table)
    {
        for (auto const & y : kvp.second)
        {
            ys.push_back(y);
        }
        yss.push_back(move(ys));
    }
    // return Polynomial(...);
}

int main(int argc, char ** argv)
{
    if (argv[1] == nullptr)
    {
        cout << "Usage: " << argv[0] << " filename" << endl;
        cout << "In file:" << '\n'
             << "    x_1 y_1,1 y_1,2 ... y_1,n_1" << '\n'
             << "    x_2 y_2,1 y_2,2 ... y_2,n_2" << '\n'
             << "       ... " << '\n'
             << "    x_m y_m,1 y_m,2 ... y_m,n_m" << '\n'
             << "where n_1 ... n_m may vary, but at least one," << '\n'
             << "x_k all distinct (equal x_k overrides previous ones)" << endl;
        return 0;
    }
    if (hermite_read(argv[1]))
        hermite_show();
    else
    {
        cerr << "Failed reading from file " << argv[1] << "!" << endl;
    }
    return 0;
}