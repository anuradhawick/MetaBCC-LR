#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

u_int64_t revComp(u_int64_t x, size_t sizeKmer=3)
{
    u_int64_t res = x;

    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*( 32 - sizeKmer))) ;
}

vector<double> getKmers(string &line, int size = 3)
{
    double total = 0;
    long  len = 0;
    u_int64_t val = 0;
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

    for (int i = 0; i < (int)line.length(); i++)
    {        
        if (!(line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }

        val = (val << 2);
        val = val & 63;
        val += (line[i] >> 1 & 3);
        len++;

        if (len == size)
        {
            // use val as the kmer for counting
            len--;  
            kmers[min(val, revComp(val, 3))]++;
            total++;          
        }
    }

    for (map<int, double>::iterator it = kmers.begin(); it != kmers.end(); ++it)
    {
        stats.push_back((it->second / max(1.0, total)));
    }

    return stats;
}

void processLinesBatch(vector<string> &linesBatch, string &outputPath, int threads)
{
    vector<vector<double>> results(linesBatch.size());
    ofstream output;
    output.open(outputPath, ios::out | ios::app);
    string o = "";

    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
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
        if (lineNum % 2 != 1)
        {
            lineNum++;
            continue;
        }
        else
        {
            batch.push_back(line);
        }
        lineNum++;
        if (batch.size() == 100000)
        {
            processLinesBatch(batch, outputPath, threads);
            batch.clear();
        }
    }

    processLinesBatch(batch, outputPath, threads);

    myfile.close();
    batch.clear();

    return 0;
}