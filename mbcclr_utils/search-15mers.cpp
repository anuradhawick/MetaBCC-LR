#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include "io_utils.h"

using namespace std;

queue<string> reads_queue;
mutex mux;
condition_variable condition;
volatile bool terminate_threads;

u_int64_t revComp(u_int64_t x, size_t sizeKmer = 15)
{
    u_int64_t res = x;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2 * (32 - sizeKmer)));
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

vector<atomic<u_int32_t>> readKmerFile(string filename)
{
    ifstream input;
    u_int64_t size;
    input.open(filename);

    input.read(reinterpret_cast<char *>(&size), sizeof(size));
    vector<atomic<u_int32_t>> kmers(size);

    input.read(reinterpret_cast<char *>(kmers.data()), kmers.size() * sizeof(kmers[0]));
    input.close();

    return kmers;
}

double *processLine(string &line, vector<atomic<u_int32_t>> &allKmers, long bin_size, int bins)
{
    double *counts = new double[bins];
    long sum = 0, count, pos, len = 0;
    u_int64_t val = 0;

    // to avoid garbage memory
    for (int i = 0; i < bins; i++)
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
            count = count < 2 ? 0: count;
            pos = (count / bin_size) - 1;

            if (count <= bin_size)
            {
                counts[0]++;
            }
            else if (pos < bins && pos > 0)
            {
                counts[pos]++;
            }
            else
            {
                counts[bins - 1]++;
            }
            sum++;
        }
    }

    if (sum > 0)
    {
        for (int i = 0; i < bins; i++)
        {
            counts[i] /= sum;
            if (counts[i] < 1e-4)
            {
                counts[i] = 0;
            }
        }
    }

    return counts;
}

void processLinesBatch(vector<string> &linesBatch, vector<atomic<u_int32_t>> &allKmers, string &output_path, int threads, long bin_size, int bins)
{
    vector<double *> batchAnswers(linesBatch.size());

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint i = 0; i < linesBatch.size(); i++)
    {
        batchAnswers[i] = processLine(linesBatch[i], allKmers, bin_size, bins);
    }

    ofstream output;
    output.open(output_path, ios::out | ios::app);
    string results = "";

    for (uint i = 0; i < linesBatch.size(); i++)
    {
        for (int j = 0; j < bins; j++)
        {
            results += to_string(batchAnswers[i][j]);

            if (j < bins - 1)
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

void off_load_process(string &output, vector<atomic<u_int32_t>> &all_kmers, int &threads, long bin_size, int bins)
{
    string seq;
    vector<string> batch;

    while (true)
    {
        {
            unique_lock<mutex> lock(mux);

            while (reads_queue.size() > 0)
            {
                seq = reads_queue.front();
                batch.push_back(seq);
                reads_queue.pop();

                if (batch.size() == 10000)
                {
                    break;
                }
            }
        }

        condition.notify_all();

        if (batch.size() > 0)
        {
            processLinesBatch(batch, all_kmers, output, threads, bin_size, bins);
            batch.clear();
        }

        {
            unique_lock<mutex> lock(mux);
            if (terminate_threads && reads_queue.size() == 0)
            {
                break;
            }
        }
    }
}

void io_thread(string &file_path)
{
    SeqReader reader(file_path);
    Seq seq;
    int count = 0;

    while (reader.get_seq(seq))
    {
        {
            unique_lock<mutex> lock(mux);
            condition.wait(lock, [] { return reads_queue.size() < 50000; });
            reads_queue.push(seq.data);
        }
        count++;
        cout << "Loaded Reads " << count << "       \r" << flush;
    }

    cout << endl;

    terminate_threads = true;
}

int main(int argc, char **argv)
{
    vector<string> batch;
    string kmers_file = argv[1];
    cout << "K-Mer file " << kmers_file << endl;

    cout << "LOADING KMERS TO RAM" << endl;
    vector<atomic<u_int32_t>> kmers = readKmerFile(kmers_file);

    cout << "FINISHED LOADING KMERS TO RAM " << endl;

    string input_path = argv[2];
    string output_path = argv[3];
    int bin_size = stoi(argv[4]);
    int bins = stoi(argv[5]);
    int threads = stoi(argv[6]);

    cout << "INPUT FILE " << input_path << endl;
    cout << "OUTPUT FILE " << output_path << endl;
    cout << "THREADS " << threads << endl;
    cout << "BIN WIDTH " << bin_size << endl;
    cout << "BINS IN HIST " << bins << endl;

    thread iot(io_thread, ref(input_path));
    thread process(off_load_process, ref(output_path), ref(kmers), ref(threads), ref(bin_size), ref(bins));

    iot.join();
    process.join();

    kmers.clear();

    cout << "COMPLETED : Output at - " << output_path << endl;

    return 0;
}
