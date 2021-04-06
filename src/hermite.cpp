#include "polynomial.hpp"

#include <cstdlib>
#include <cctype> // isspace

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

typedef Polynomial<double>  Poly;

map<double, vector<double>> Table; // the table storing the input.
Poly                        P; // The generated Polynomial
vector<double>              A; // The generated pre-coefficients

bool                        Verbose = false; // Comment on the actions.
bool                        Expanded = false, Decomposite = true; // see options.
int                         Width = 20; // Width of number formatting in tabualar output.

//const char ERROR[] = "\x1b[31;1mError\x1b[0m: ";
const char ERROR[] = "Error: ";

const char INTRO[] =
    "The Interactive Hermite Interpolation Generator 0.1\n"
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
    " 'o':  Options: Parameter option and parameters..\n"
    "       Example usage: 'o v +', 'o w 20'\n"
    "       Available options:\n"
    "       'v': Verbose.\n"
    "         + : Set. Display all messages. (default after CLA processing)\n"
    "         - : Unset. Only display outputs and error messages. (default while\n"
    "             CLA processing)\n"
    "       'w': Minimum width for numbers in tabular output.\n"
    "       'p': Print polynomial decomposite, expanded, or both.\n"
    "         + : expanded, i.e. form of an x^n + ... + a_0 (default)\n"
    "         - : decomposite, i.e. form of b0 + b1(x - x1) + ... +\n"
    "                                       bn(x - x1)^k1...(x - xm)^km\n"
    "         0 : both, i.e. '<decomposite>  =  <expanded>'\n"
    " 'l':  Load. Parameter filename. Load a file. This must be done first.\n"
    "       The command exists to reload the file or view another one\n"
    "       without exiting.\n"
    "       Any following commands expect a file to be loaded.\n"
    " 'v':  Values. Show the loaded values (data that has been loaded).\n"
    " 'p':  Polynomial. Show the definig term of the interpolation polynomial.\n"
    " 'e':  Evaluate: Parameter x. Evaluates the polynomial function at x.\n"
    " 'd':  Derivative: Replaces the active polynomial by its derivative.\n"
    " 't':  Value table: Parameters l, n, u. Evaluates the polynmial between\n"
    "       l and u with n steps (i.e. n + 1 places) between\n"
    "Commands and options are indeed not really checked for existance.\n"
    "Any word starting with one of the single letters suit.\n"
    "The empty word is ignored.\n"
    "Commands and paramerets are split with any white space:\n"
    "'load data.txt  eval 0' is like\n"
    "'l data.txt e 0' or\n"
    "'lunch data.txt eat 0'\n"
    ;

bool hermite_load(char const * const filename)
{
    ifstream file(filename);
    if (!file.good())
    {
        cerr << ERROR << "Cannot open file " << filename << "!\n";
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
            cerr << ERROR << "Cannot read x value in line " << i << "!\n";
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
            cerr << ERROR << "Cannot read y value no. " << j << " in line " << i << "!\n";
            goto fail;
        }
        if (j < 2)
        {
            cerr << ERROR << "Value " << x << "in line " << i << " needs at least one value!";
            goto fail;
        }
        Table[x] = move(ys);
    }

    if (Verbose) cout << "File '" << filename << "' successfully loaded." << endl;
    return true;

  fail:
    Table.clear();
    A.clear();
    P = 0;
    cerr << "Failing load reset state." << endl;
    return false;
}

inline
void print_v(vector<double> const & v)
{
    cout << '[';
    if (!v.empty())
    {
        cout << setw(Width) << v.front();
        for (std::size_t i = 1; i < v.size(); ++i)
            cout << ',' << setw(Width) << v[i];
    }
    cout << ']';
}

void print_v(vector<vector<double>> const & v)
{
    cout << "{\n";
    if (!v.empty())
    {
        print_v(v.front());
        for (std::size_t i = 1; i < v.size(); ++i)
        {
            cout << '\n';
            print_v(v[i]);
        }
    }
    cout << "\n}";
}

vector<double> hermite_ipol()
{
    auto fac = [](auto n)
        {
            double r = 1.0;
            for (decltype(n) i = 2; i < n; ++i) r *= i;
            return r;
        };

    std::size_t n = 0; // size of the system aka. total # of ys in Table.
    vector<double> xs;
    vector<vector<double>> yss;
    for (auto const & kvp : Table)
    {
        int i = n;
        n += kvp.second.size();
        xs.resize(n, kvp.first);
        yss.resize(n, vector<double>(1, kvp.second.front()));
        for (std::size_t k = 1; k < kvp.second.size(); ++k)
        for (std::size_t j = k; j < kvp.second.size(); ++j)
            yss[i + j].push_back(kvp.second[k] / fac(k+1));
    }

    for (std::size_t i = 1; i < n; ++i)
    {
        auto l = xs.begin();
        auto r = l + i;
        for (int j = i; r != xs.end(); ++l, ++r, ++j)
        {
            if (*r == *l) continue;
            double y = (yss[j][i-1] - yss[j-1][i-1]) / (*r - *l);
            yss[j].push_back(y);
        }
    }
    if (Verbose) { print_v(yss); cout << "\nCoefficients successfully calulated." << endl; }

    A = vector<double>(n);
    for (std::size_t i = 0; i < n; ++i) A[i] = yss[i][i];
    return xs;
}

bool interactive(istream & in)
{
    auto ex_table = []() -> bool
        {
            if (Table.empty()) cerr << ERROR << "No file has been loaded." << endl;
            return !Table.empty();
        };

    string command, opt;
    if (in >> command)
    switch (command[0])
    {
      case 'h': cout << HELP << endl; break;
      case 'q': if (Verbose) cout << "Leaving." << endl; exit(0);
      case 'o':
        if (in >> opt)
        {
            switch (opt[0])
            {
              case 'v': if (in >> opt)
                {
                    if (opt[0] == '+') { Verbose = true; cout << "Verbose on." << endl; }
                    else if (opt[0] == '-') Verbose = false;
                    else
                    {
                        cerr << ERROR << "Illegal option parameter; only '+' or '-' allowed.\n";
                        return false;
                    }
                }
                else
                {
                    cerr << ERROR << "Failed to read option parameter.\n";
                    in.setstate(ios::goodbit);
                    return false;
                }
                break;
              case 'a': if (in >> opt)
                {
                    if (opt[0] == '0' || opt[0] == '+' || opt[0] == '-')
                    {
                        Expanded = (opt[0] == '0' || opt[0] == '+');
                        Decomposite = (opt[0] == '0' || opt[0] == '-');
                    }
                    else
                    {
                        cerr << ERROR << "Illegal option parameter; only '+' or '-' allowed.\n";
                        return false;
                    }
                }
                else
                {
                    cerr << ERROR << "Failed to read option parameter.\n";
                    in.setstate(ios::goodbit);
                    return false;
                }
                break;
              case 'w':
                if (in >> Width)
                {
                    if (Width < 1)
                    {
                        Width = 20;
                        if (Verbose) cout << "Negative number. Width set to 20 (default)." << endl;
                    }
                    else
                        if (Verbose) cout << "Negative number. Width set to " << Width << "." << endl;
                }
                else
                {
                    if      (in.eof())  cerr << ERROR << "No option parameter.\n";
                    else if (in.fail()) cerr << ERROR << "Illegal option parameter format. Must be positive number.\n";
                    in.setstate(ios::goodbit);
                    return false;
                }
                break;
              default:
                cerr << ERROR << "Unknown option '" << opt << "'.\n";
                return false;
            }
        }
        else
        {
            cerr << ERROR << "Failed to read option.\n";
            in.setstate(ios::goodbit);
            return false;
        }
        break;
      case 'v': if (ex_table())
        {
            for (auto const & kvp : Table)
            {
                cout << setw(Width) << kvp.first << " : ";
                for (auto const & y : kvp.second) cout << setw(Width) << y;
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
                auto x = hermite_ipol();
                if (!A.size())
                {
                    cerr << ERROR << "No coefficents generated. Sorry, this is an implementation error and not your fault." << endl;
                    return false;
                }
                Poly N = 1.0;
                auto n = [&N, &x](int k) { return N *= Poly(1, 1) - x[k]; };

                P = A.front();
                for (std::size_t i = 1; i < A.size(); ++i) P += A[i] * n(i - 1);
                if (Verbose) cout << "Polynomial successfully calculated." << endl;
            }
            else
            {
                cerr << ERROR << "Failed to read filename.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
        break;
      case 'p': if (ex_table())
        {
            //TODO: Use A to print P in the form of a0 + a1(x - x0) + a2(x - x0)(x - x1) + ... + an(x - x1)...(x - xn)
            if (Expanded)
            {
                cout << A[0];
                auto itT = Table.begin();
                auto itA = A.begin();
                cout << *itA;
                vector<pair<Poly, std::size_t>> Qds = { { Poly(1, 1) - itT->first, 1 } };
                auto inc = [&itT, &Qds]()
                    {
                        if (Qds.back().second < itT->second.size())
                            ++Qds.back().second;
                        else
                            Qds.push_back(pair{ Poly(1, 1) - (++itT)->first, 1 });
                    };
                for (++itA; itA != A.end(); inc(), ++itA)
                {
                    if (*itA == 0.0) continue;
                    if (*itA <  0.0) cout << " - " << -*itA;
                    else             cout << " + " <<  *itA;
                    for (auto const & Qd : Qds) // Q: Poly, d: degree
                    {
                        auto const & Q = Qd.first;
                        auto const & d = Qd.second;
                        if (Q.is_monomial()) cout << ' ' <<        Q;
                        else                 cout << ' ' << '(' << Q << ')';
                        if (d > 1) cout << '^' << d;
                    }
                }
                if (Decomposite) cout << "  =  ";
                else             cout << '\n';
            }
            if (Decomposite) cout << P << endl;
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
                cerr << ERROR << "Failed to interpret the parameter as IEEE-754 double number.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
        break;
      case 'd': if (ex_table())
        {
            // Adapt A.
            P = P.deriv();
            cout << "Polynmial replaced by its derivative." << endl;
        }
        break;
      case 't': if (ex_table())
        {
            double l, u;
            int n;
            if ((in >> l >> n >> u) && n > 0)
            {
                double const d = (u - l) / n;
                for (int i = 0; i <= n; ++i)
                {
                    double const x = l + i * d;
                    cout << setw(Width) << x << setw(Width) << P(x) << '\n';
                }
                cout << endl;
            }
            else
            {
                cerr << ERROR << "Failed to read Paramters. Must be double, positive integer, double.\n";
                in.setstate(ios::goodbit);
                return false;
            }
        }
        break;
      default:
        cerr << ERROR << "Unknown command.\n";
    }
    return !in.eof();
}

int main(int argc, char ** argv)
{
    if (argc < 2) cout << INTRO << endl;

    stringstream str;
    while (--argc && ++argv) str << *argv << ' ';
    while (interactive(str)) { }

    if (!str.eof())
    {
        cerr << "Some error occoured. Leaving." << endl;
        return 1;
    }

    Verbose = true;
    str = stringstream(); // reinitialize

    string line;
    do
    {
        cout << "hermite> " << flush;
        if (!getline(cin, line)) break;
        while (!line.empty() && isspace(line.back())) line.pop_back();
    }
    while (line.empty() || interactive(reinterpret_cast<stringstream &>(str << line << ' ')));

    return 0;
}