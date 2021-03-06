#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

u_int64_t revComp(u_int64_t x, size_t sizeKmer=15)
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

vector<long> splitLine(string &line)
{
    vector<long> vec(2);
    string tmp = "";

    for (int i = 0; i < (int)line.length(); i++)
    {
        if (line[i] == ',' || line[i] == ' ')
        {
            vec[0] = stol(tmp);
            tmp = "";
            continue;
        }
        else
        {
            tmp += line[i];
        }
    }

    vec[1] = stol(tmp);

    return vec;
}

long readKmerFile(string filename, vector<u_int32_t> &kmers)
{
    ifstream myfile(filename);
    string line;
    vector<long> v;
    long count = 0;
    while (getline(myfile, line))
    {
        v = splitLine(line);
        kmers[v[0]] = v[1];
        kmers[revComp(v[0])] = v[1];
        count ++;
    }

    myfile.close();

    return count;
}

double *processLine(string &line, vector<u_int32_t> &allKmers, long bin_size = 10)
{
    double *counts = new double[32];
    long sum = 0, count, pos, len = 0;
    u_int64_t val = 0;

    // to avoid garbage memory
    for (int i = 0; i < 32; i++)
    {
        counts[i] = 0;
    }

    for (size_t i = 0; i < line.length(); i++)
    {
        if (!(line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }

        val = (val << 2);
        val = val & 1073741823;
        val += (line[i] >> 1 & 3);
        len++;

        if (len == 15)
        {
            // use val as the kmer for counting
            len--;  
            count = allKmers[(long)val];
            pos = (count / bin_size) - 1;

            if (count <= bin_size)
            {
                counts[0]++;
                
            }
            else if (pos < 32 && pos > 0)
            {
                counts[pos]++;
            }
            sum++;          
        }
    }

    for (int i = 0; i < 32; i++)
    {
        counts[i] /= sum;
        if (counts[i] <  1e-4)
        {
            counts[i] = 0;
        }
    }

    return counts;
}

void processLinesBatch(vector<string> &linesBatch, vector<u_int32_t> &allKmers, string &outputPath, int threads, long bin_size = 10)
{
    vector<double *> batchAnswers(linesBatch.size());

    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint i = 0; i < linesBatch.size(); i++)
    {
        batchAnswers[i] = processLine(linesBatch[i], allKmers, bin_size);
    }

    ofstream output;
    output.open(outputPath, ios::out | ios::app);
    string results = "";

    for (uint i = 0; i < linesBatch.size(); i++)
    {
        for (int j = 0; j < 32; j++)
        {
            results += to_string(batchAnswers[i][j]);

            if (j < 32 - 1)
            {
                results += " ";
            }
        }

        results += "\n";

        // releasing pointer memory
        delete[] batchAnswers[i];
    }

    output << results;
    batchAnswers.clear();
    output.close();
}

int main(int argc, char ** argv)
{
    vector<u_int32_t> kmers(1073741824, 0);    
    vector<string> batch;
    long lineNum = 0, count;

    string kmersFile =  argv[1];
    cout << "K-Mer file " << kmersFile << endl;

    cout << "LOADING KMERS TO RAM" << endl;
    count = readKmerFile(kmersFile, kmers);

    cout << "FINISHED LOADING KMERS TO RAM " << count << endl;

    string inputPath =  argv[2];
    string outputPath =  argv[3];
    int bin_size =  max(stoi(argv[4]), 10);
    int threads = 8;

    if (argv[5] != NULL) 
        threads = stoi(argv[5]);

    cout << "INPUT FILE " << inputPath << endl;
    cout << "OUTPUT FILE " << outputPath << endl;
    cout << "THREADS " << threads << endl;
    cout << "BIN WIDTH " << bin_size << endl;

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
            processLinesBatch(batch, kmers, outputPath, threads, bin_size);
            batch.clear();
        }
    }
    processLinesBatch(batch, kmers, outputPath, threads);

    myfile.close();
    batch.clear();

    kmers.clear();

    cout << "COMPLETED : Output at - " << outputPath << endl;

    return 0;
}
