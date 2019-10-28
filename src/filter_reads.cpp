#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char ** argv)
{
    string line1, line2, line3, line4, truth, 
    readsPath = argv[1],
    outputPath = argv[2],
    type = argv[2];
    
    ifstream infile(readsPath);
    ofstream outfile(outputPath, ios::out);

    // TODO update to handle multi line fasta
    
    while (getline(infile, line1))
    {
        getline(infile, line2);

        if (type == "fq") {
            getline(infile, line3);
            getline(infile, line4);
        }
        
        if (line2.length() >= 1000)
        {
            if (type == "fq") {
                outfile << "<" << line1.substr(1, line1.length()) << endl;
            } else {
                outfile << line1 << endl;
            }
            outfile << line2 << endl;
        }
    }

    outfile.close();
    infile.close();

    return 0;
}