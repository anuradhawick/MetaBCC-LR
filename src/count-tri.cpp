#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <regex>

using namespace std;


int val(char &c)
{
    if (c == 'A')
    {
        return 0;
    }
    else if (c == 'C')
    {
        return 1;
    }
    else if (c == 'T')
    {
        return 2;
    }
    else
    {
        return 3;
    }
}

string revCmpl(string kmer)
{
    string reverseComplement = "";
    static map<char, char> cmpl = {
        {'A', 'T'},
        {'G', 'C'},
        {'C', 'G'},
        {'T', 'A'}};

    for (int i = kmer.length() - 1; i >= 0; i--)
    {
        reverseComplement += cmpl[kmer[i]];
    }

    return reverseComplement;
}

long toDeci(string str, int base = 4)
{
    int llen = str.length(), power = 1;
    long num = 0;
    for (int i = llen - 1; i >= 0; i--)
    {
        if (val(str[i]) >= base)
        {
            cout << "Invalid Number";
            return -1;
        }
        num += val(str[i]) * power;
        power = power * base;
    }

    return num;
}

vector<double> getKmers(string &readLine, int size = 3)
{
    static regex validKmers("^[CAGT]+$");
    double total = 0;
    map<int, double> kmers = {
        {0, 0},
        {1, 0},
        {3, 0},
        {2, 0},
        {4, 0},
        {5, 0},
        {7, 0},
        {6, 0},
        {12, 0},
        {13, 0},
        {15, 0},
        {8, 0},
        {9, 0},
        {11, 0},
        {16, 0},
        {17, 0},
        {19, 0},
        {20, 0},
        {21, 0},
        {23, 0},
        {28, 0},
        {29, 0},
        {24, 0},
        {25, 0},
        {41, 0},
        {49, 0},
        {45, 0},
        {53, 0},
        {37, 0},
        {33, 0},
        {32, 0},
        {36, 0}
    };
    vector<double> stats;

    for (int i = 0; i < (int)readLine.length() - size - 1; i++)
    {
        //ignore kmers with non ACGT characters
        if (!regex_match(readLine.substr(i, size),  validKmers)) {
            continue;
        }
        kmers[min(toDeci(readLine.substr(i, size)), toDeci(revCmpl(readLine.substr(i, size))))]++;
        total++;
    }

    for (map<int, double>::iterator it = kmers.begin(); it != kmers.end(); ++it)
    {
        stats.push_back((it->second / max(1.0, total)));
    }

    return stats;
}

void processLinesBatch(vector<string> &linesBatch, string &outputPath)
{
    vector<vector<double>> results(linesBatch.size());
    ofstream output;
    output.open(outputPath, ios::out | ios::app);
    string o = "";

    #pragma omp parallel for num_threads(8)
    for (size_t i = 0; i < linesBatch.size(); i++)
    {
        results[i] = getKmers(linesBatch[i]);
    
    }

    for (size_t i = 0; i < linesBatch.size(); i++)
    {
        for (double d: results[i])
        {
            o += to_string(d);
            o += " ";
        }
        o += "\n";
    }

    output << o;
    results.clear();
    output.close();
}

int main(int argc, char ** argv)
{
    vector<string> batch;
    long lineNum = 0;

    string inputPath = argv[1];
    string outputPath = argv[2];
    int threads = 8;
    if (argv[3]!= NULL) 
        threads = stoi(argv[3]);

    cout << "INPUT FILE " << inputPath << endl;
    cout << "OUTPUT FILE " << outputPath << endl;
    cout << "THREADS " << threads << endl;

    ifstream myfile(inputPath);
    string line;

    ofstream output;
    output.open(outputPath, ios::out);
    output.close();

    while (getline(myfile, line))
    {
        if (lineNum % 4 != 1)
        {
            lineNum++;
            continue;
        }
        else
        {
            batch.push_back(line);
        }
        lineNum++;
        if (batch.size() == 10000)
        {
            processLinesBatch(batch, outputPath);
            batch.clear();
        }
    }

    processLinesBatch(batch, outputPath);

    myfile.close();
    batch.clear();

    return 0;
}