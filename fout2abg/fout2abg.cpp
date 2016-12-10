/*
 * Copyright (c) 2016 Hailin Su, ISU NBCEC
 *
 * This file is part of SNPipeline.
 *
 * SNPipeline is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SNPipeline is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Lesser Public License
 * along with SNPipeline.  If not, see <http://www.gnu.org/licenses/>.
 */

//
//  fout2abg.cpp
//  libHailins
//
//  Created by Hailin SU on 3/31/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

//  fout2abg

#include "../appInfo.h"
#include <vector>

using namespace std;
namespace libcbk
{
    class fout2abg : public appInfo
    {
    public:
        fout2abg(int argc, const char * argv[]);
        vector<string> vecCoding;
    private:
        void theMain(int argc, const char * argv[]);
    };

    fout2abg::fout2abg(int argc, const char * argv[])
    {
        progName    = "fout2abg";
        version     = "2014-04-18 11:10 CDT";
        created     = "2014-03-31";
        updated     = "2014-04-18";
        descrip     = "Converts FImpute output file to ab-genotype file.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-03-31 init
         2014-04-18 checks for the consistance of map file and genotype file
         */
    }

    void fout2abg::theMain(int argc, const char * argv[])
    {
        makeSure(argc == 5, argv[0], "snp_info.txt genotypes_imp.txt -o ab-genotype.output.abg");
        makeSure(argv[3] == string("-o"), argv[0], "snp_info.txt genotypes_imp.txt -o ab-genotype.output.abg");

        I("s", 1, "Reading SNP_IDs...");

        ifstream infs (argv[1]);
        if(!infs)
        {
            E("ss", "File not found:", argv[1]);
            exit(-1);
        }
        ofstream ofs(argv[4]);
        if(!ofs)
        {
            infs.close();
            E("ss", "Could not open file to write:", argv[4]);
            exit(-3);
        }
        string stmp, aLine("#SNP_IDs");
        unsigned long long fsp1(0.01*tellFileSize(argv[1]));
        if(fsp1<1)
            fsp1=1;
        getline(infs, stmp);    // get rid of header line
        unsigned linesR(0);
        while(infs>>stmp)   // get a SNP_ID
        {
            aLine+="#";
            aLine+=stmp;
            getline(infs, stmp);    // get rid of the rest of the line
            if(++linesR%200==0)
                progShow('L', infs.tellg()/fsp1);
        }
        infs.close();
        const unsigned numberSNPs (linesR);
        // write out the header
        ofs << aLine << endl;
        progClear();

        I("s", 1, "Converting genotypes_imp.txt to abg file...");
        infs.open(argv[2]);
        if(!infs)
        {
            ofs.close();
            E("ss", "File not found:", argv[2]);
            exit(-2);
        }
        fsp1=0.01*tellFileSize(argv[2]);
        if(fsp1<1)
            fsp1=1;
        linesR=0;
        getline(infs, stmp);    // get rid of header line
        unsigned i;
        vecCoding.resize(6);
        vecCoding[0]="AA";
        vecCoding[2]="BB";
        vecCoding[1]=vecCoding[3]=vecCoding[4]="AB";
        vecCoding[5]="--";
        while(infs>>aLine)   // read the animalID
        {
            infs>>stmp>>stmp;     // the chip column and genotype column
            if(stmp.length()!=numberSNPs)
            {
                E("si", "numberSNPs read from snp_info.txt :", numberSNPs);
                E("si", "numberSNPs read from genotype     :", stmp.length());
                E("ss", "in line starting with             :", aLine.c_str());
                exit(-2);
            }
            for(i=0; i<numberSNPs; ++i)
            {
                aLine+=" ";
                aLine+=vecCoding[stmp[i]-'0'];
            }
            ofs << aLine << endl;
            if(++linesR%200==0)
                progShow('C', infs.tellg()/fsp1);
        }
        progClear();
        I("ss", 1, argv[4], "generated.");
        infs.close();
        ofs.close();
    }
}

int main(int argc, const char * argv[])
{
    libcbk::fout2abg(argc, argv);
    return 0;
}