#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char ** argv)
{
    string line1, line2, line3, line4, truth, 
    readsPath = argv[1],
    outputPath = argv[2],
    type = argv[3],
    truthPath, truthOutputPath;
    ofstream outfileTruth;
    ifstream truthfile;

    if (argc == 6)
    {
        truthPath = argv[4];
        truthfile.open(truthPath);
        truthOutputPath = argv[5];
        outfileTruth.open(truthOutputPath, ios::out);
    }
    
    ifstream infile(readsPath);
    ofstream outfile(outputPath, ios::out);

    // TODO update to handle multi line fasta
    
    while (getline(infile, line1))
    {
        getline(infile, line2);

        if (argc == 6)
        {
            getline(truthfile, truth);
        }

        if (type == "fq") {
            getline(infile, line3);
            getline(infile, line4);
        }
        
        if (line2.length() >= 1000)
        {
            outfile << ">" << line1.substr(1, line1.length()) << endl;
            outfile << line2 << endl;
        
            if (argc == 6)
            {
                outfileTruth << truth << endl;
            }
        }
    }

    outfile.close();
    infile.close();

    if (argc == 6)
    {
        outfileTruth.close();
        truthfile.close();
    }

    return 0;
}