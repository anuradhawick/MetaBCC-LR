#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char ** argv)
{
    string line1, line2, line3, line4, truth, 
    readsPath = argv[1],
    outputPath = argv[2];
    
    ifstream infile(readsPath);
    ofstream outfile(outputPath, ios::out);
    
    while (getline(infile, line1))
    {
        getline(infile, line2);
        getline(infile, line3);
        getline(infile, line4);

        if (line2.length() >= 1000)
        {
            outfile << line1 << endl;
            outfile << line2 << endl;
            outfile << line3 << endl;
            outfile << line4 << endl;
        }
    }

    outfile.close();
    infile.close();

    return 0;
}