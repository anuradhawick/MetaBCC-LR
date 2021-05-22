#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
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

void writeKmerFile(string filename, vector<atomic<u_int32_t>> &kmers)
{
    ofstream output;
    output.open(filename, ios::out);
    u_int64_t size = kmers.size();
    output.write(reinterpret_cast<char const*>(&size), sizeof(size));
    output.write(reinterpret_cast<char const*>(kmers.data()), kmers.size() * sizeof(kmers[0]));
    output.close();
}

vector<atomic<u_int32_t>> readKmerFile(string filename)
{
    ifstream input;
    u_int64_t size;
    input.open(filename);

    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    vector<atomic<u_int32_t>> kmers(size);

    input.read(reinterpret_cast<char*>(kmers.data()), kmers.size() * sizeof(kmers[0]));
    input.close();

    return kmers;
}

void processLine(string &line, vector<atomic<u_int32_t>> &all_kmers)
{
    long len = 0;
    u_int64_t val = 0, rev = 0;
    u_int32_t oval;

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
            // record original value from index
            oval = all_kmers[val];
            
            // CAS
            while (!all_kmers[val].compare_exchange_weak(oval, oval + 1))
            {
            };

            // reverse complement recording
            rev = revComp(val);
            oval = all_kmers[rev];
            
            // CAS
            while (!all_kmers[rev].compare_exchange_weak(oval, oval + 1))
            {
            };
        }
    }
}

void processLinesBatch(vector<string> &batch, vector<atomic<u_int32_t>> &all_kmers, int threads)
{
    vector<vector<u_int64_t>> batch_results(batch.size());

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (size_t i = 0; i < batch.size(); i++)
    {
        processLine(batch[i], all_kmers);
    }

    batch_results.clear();
}

void off_load_process(string &output, vector<atomic<u_int32_t>> &all_kmers, int &threads)
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

                if (batch.size() == 1000)
                {
                    break;
                }
            }
        }

        condition.notify_all();

        if (batch.size() > 0)
        {
            processLinesBatch(batch, all_kmers, threads);
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
            condition.wait(lock, [] { return reads_queue.size() < 10000; });
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
    vector<atomic<u_int32_t>> kmers(1073741824);

    string input_path = argv[1];
    string output_path = argv[2];
    int threads = stoi(argv[3]);

    cout << "INPUT FILE " << input_path << endl;
    cout << "OUTPUT FILE " << output_path << endl;
    cout << "THREADS " << threads << endl;

    thread iot(io_thread, ref(input_path));
    thread process(off_load_process, ref(output_path), ref(kmers), ref(threads));

    iot.join();
    process.join();

    cout << "WRITING TO FILE" << endl;

    writeKmerFile(output_path, kmers);
    
    kmers.clear();
    cout << "COMPLETED : Output at - " << output_path << endl;

    return 0;
}
