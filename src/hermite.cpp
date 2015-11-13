// Compile with clang++-3.0 -std=c++11

#include "polynomial.hpp"

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <utility> // pair
#include <map>
#include <vector>

using namespace std;
using namespace Modulus;

using Poly = Polynomial<double>;

map<double, vector<double>> Table;
bool                        Verbose = false;
Poly                        P;

const char INTRO[] =
    "Hermite Interpolation Generator 0.1\n"
    "by Q. F. Schroll @ github.com/Bolpat\n"
    "NO WARRENTY, USE FOR ANY PURPOSE\n\n"
    
    "Type in 'help' or 'h' for further information and how to use."
    ;

const char HELP[] =
    "This program is interactive, i.e. it waits for you to type commands in.\n"
    "Any command-line arguments will be trated as they were typed commands\n"
    "in the interactive mode (only lacks some less relevant verbose output).\n\n"
    
    "Before any form of evaluation, you have to load a file.\n"
    "That file needs a form like:\n"
    "    x_1   y_1,1 y_1,2 ... y_1,n_1\n"
    "    x_2   y_2,1 y_2,2 ... y_2,n_2\n"
    "       ... \n"
    "    x_m   y_m,1 y_m,2 ... y_m,n_m\n"
    "where n_1, ..., n_m may vary, but all at least one,\n"
    "and x_k preferably distinct (equal x_k override previous ones)\n"
    "# indicates line comments if it is the first character in that line.\n"
    "Empty lines are ok but ignored.\n\n"
    
    "List of interactive commands:\n"
    " 'h':  Help. Shows this.\n"
    " 'q':  Quit. Immediately exit the program.\n"
    " 'l':  Load. Parameter filename. Load a file. This must be done first.\n"
    "       The command exists to reload the file or view another one\n"
    "       without exiting.\n"
    "       Any following commands expect a file to be loaded.\n"
    " 'v':  Values. Show the loaded values (data that has been loaded).\n"
    " 'p':  Polynomial. Show the definig term of the interpolation polynomial.\n"
    " 'e':  Evaluate: Parameter x. Evaluates the polynomial function at x.\n"
    " 'd':  Derivative: Replaces the active polynomial by its derivative.\n"
    " 't':  Value table: Parameters l, n, u. Evaluates the polynmial between l and u\n"
    "       with n + 1 steps between\n"
    "Commands are indeed not really checked for existance. Any word starting with one\n"
    "of the single letters suit. The empty word is ignored.\n"
    "Just split commands and paramerets:\n"
    "'load data.txt  eval 0' is like\n"
    "'l data.txt e 0' or\n"
    "'lunch data.txt eat 0'\n"
    ;

bool hermite_load(char const * const filename)
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
        if (line == "" || line[0] == '#') continue;
        istringstream iss(line);
        if (!(iss >> x))
        {
            cerr << "Error: Cannot read x value in line " << i << "!\n";
            goto fail;
        }
        int j = 1;
        while (iss >> y)
        {
            ys.push_back(y);
            ++j; // expect j-th value next.
        }
        if (!iss.eof())
        {
            cerr << "Error: Cannot read y value no. " << j << " in line " << i << "!\n";
            goto fail;
        }
        if (j < 2)
        {
            cerr << "Error: " << x << "in line " << i << " needs at least one value!";
            goto fail;
        }
        Table[x] = move(ys);
    }
    
    if (Verbose) cout << "File '" << filename << "' successfully loaded." << endl;
    return true;
    
  fail:
    Table.clear();
    if (Verbose) cout << "Failing load resetted state." << endl;
    return false;
}

inline
void print_v(vector<double> const & v)
{
    cout << '[';
    if (!v.empty())
    {
        cout << setw(3) << v.front();
        for (int i = 1; i < v.size(); ++i)
            cout << ',' << setw(3) << v[i];
    }
    cout << ']';
}

void print_v(vector<vector<double>> const & v)
{
    cout << "{\n";
    if (!v.empty())
    {
        print_v(v.front());
        for (int i = 1; i < v.size(); ++i)
        {
            cout << '\n';
            print_v(v[i]);
        }
    }
    cout << "\n}";
}

int fac(int n)
{
    int r = 1;
    for (int i = 2; i < n; ++i) r *= i;
    return r;
}

pair<vector<double>,vector<double>> hermite_ipol()
{
    int n = 0; // size of the system aka. total # of ys in Table.
    vector<double> xs;
    vector<vector<double>> yss;
    for (auto const & kvp : Table)
    {
        int i = n;
        n += kvp.second.size();
        xs.resize(n, kvp.first);
        yss.resize(n, vector<double>(1, kvp.second.front()));
        for (int k = 1; k < kvp.second.size(); ++k)
        for (int j = k; j < kvp.second.size(); ++j)
            yss[i + j].push_back(kvp.second[k] / fac(k+1));
    }
    
    if (Verbose) { print_v(yss); cout << '\n'; }
    for (int i = 1; i < n; ++i)
    {
        auto l = xs.begin();
        auto r = l + i;
        for (int j = i; r != xs.end(); ++l, ++r, ++j)
        {
            if (*r == *l) continue;
            double y = (yss[j][i-1] - yss[j-1][i-1]) / (*r - *l);
            yss[j].push_back(y);
            if (Verbose) { print_v(yss); cout << '\n'; }
        }
    }
    
    vector<double> as(n);
    
    for (int i = 0; i < n; ++i)
        as[i] = yss[i][i];
    
    if (Verbose) cout << "Coefficients successfully calulated." << endl;
    
    return make_pair(move(xs), move(as));
}

bool ex_table()
{
    if (Table.empty()) cerr << "Error: No file has been loaded." << endl;
    return !Table.empty();
}

bool interactive(istream & in)
{
    string command;
    if (in >> command) switch (command[0])
    {
        case 'h': cout << HELP << endl; break;
        case 'q': if (Verbose) cout << "Leaving." << endl; exit(0);
        case 'v': if (ex_table())
        {
            for (auto const & kvp : Table)
            {
                cout << setw(8) << kvp.first << " : ";
                for (auto const & y : kvp.second) cout << setw(8) << y;
                cout << '\n';
            }
            cout << endl;
        }
            break;
        case 'l':
        {
            string s;
            if ((in >> s) && hermite_load(s.c_str()))
            {
                auto x_a       = hermite_ipol();
                auto const & x = x_a.first;
                auto const & a = x_a.second;
                if (a.size() == 0)
                {
                    cerr << "Error: no coefficents generated. This is an internal error and not your fault." << endl;
                    return false;
                }
                Poly N = 1;
                P = a.front() * N;
                for (int i = 1; i < a.size(); ++i) P += a[i] * (N *= Poly(1, 1) - x[i-1]);
                if (Verbose) cout << "Polynomial successfully calculated." << endl;
            }
            else
            {
                cerr << "Error: Failed to read filename.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
            break;

        case 'p': if (ex_table())
        {
            cout << P << endl;
        }
            break;
        case 'e': if (ex_table())
        {
            double x;
            if (in >> x)
            {
                cout << P(x) << endl;
            }
            else
            {
                cerr << "Error: Failed to interpret the parameter as IEEE-754 double number.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
            break;
        case 'd': if (ex_table())
        {
            P = P.deriv();
            if (Verbose) cout << "Polynmial replaced by its derivative." << endl;
        }
            break;
        case 't': if (ex_table())
        {
            double l, u;
            int s;
            if ((in >> l >> s >> u) && s > 0)
            {
                double d = (u - l) / s;
                for (double x = l; x <= u; x += d)
                    cout << setw(20) << x << setw(20) << P(x) << '\n';
                cout << endl;
            }
            else
            {
                cerr << "Error: Failed to read Paramters. Must be double, positive integer, double.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
            break;
    };
    return !in.eof();
}

int main(int argc, char ** argv)
{
    if (argc < 2) cout << INTRO << endl;
    
    stringstream str;
    
    while (--argc && ++argv) str << *argv << ' ';
    while (interactive(str))
        ;

    if (!str.eof())
    {
        cerr << "Some error occoured. Leaving." << endl;
        return 1;
    }
    
    Verbose = true;
    
    do cout << "> " << flush; while (interactive(cin));
    
    return 0;
}