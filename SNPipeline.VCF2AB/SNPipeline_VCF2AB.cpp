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
//  SNPipeline_VCF2AB.cpp
//  SNPipeline
//
//  Created by Hailin SU on 2/26/14.
//

#include "SNPipeline_VCF2AB.h"
#include <iomanip>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)
#define L(...)              if(vb) lg.generalLog("I", __VA_ARGS__)

using namespace std;

namespace libcbk
{
    void SNPipeline_VCF2AB::usage(const string& argv0)
    {
        cout << endl
        << "     version: " << version << endl
        << "####################################################################" << endl
        << "#######" << endl
        << "###                 " << progName << endl
        << "####" << endl
        << "############" << endl
        << "##                              by hailins@iastate.edu" << endl
        << "#####                           created " << created << endl
        << "####                            updated " << updated << endl
        << "#######" << endl
        << "####################################################################" << endl
        << endl << "    " << descrip << endl
        << endl
        << "usage:" << endl
        << argv0 << " VCF4.1_File" << endl << endl
        << endl     << endl;
    }

    SNPipeline_VCF2AB::SNPipeline_VCF2AB(const int argc, const char * argv[])
    {
        progName    = "SNPipeline.VCF2AB";
        version     = "2014-02-26 11:33 CDT";
        created     = "2014-02-26";
        updated     = "N/A";
        descrip     = "SNPipeline - setp ?. Convert phased/unphased, but must be imputed VCF4.1 file to merged/unmerged AB genotype format.";
        //        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-02-11 11:33 AM init
         */
    }

    void SNPipeline_VCF2AB::theMain(const int argc, const char * argv[])
    {
        // **************** parsing parameters
        if(argc!=2)
        {
            usage(argv[0]);
            exit (-1);
        }
        //string vcfFileName (argv[1]);
        VCF_File v(argv[1]);

        if(!v.isValid)
            return;

        //                           0     1     2     3
        //        const string vecUC2STR[4] {"--", "AA", "BB", "AB"};
        const string vecUC2STR[] = {"--", "AA", "BB", "AB"};

        const string outFile ("ab-genotype.from." + v.fileName + ".abg");

        ofstream ofs (outFile.c_str());

        if(!ofs)
        {
            E("s", "ERROR: cannot open file to write into.");
            exit (-30);
        }

        // TODO: write out the ab-genotype file, from SNPipeline/data.cpp
        unsigned long long i, j;
        const unsigned long numberSamples (v.vecSampleIDs.size());
        const unsigned long numberSNPs (v.vecSNPNames.size());
        bool vb(true);

        L("s", 1, "Start unpacking genotypes and writing them out...");

        // write out the SNP names for the whole file
        ofs << "#SNP_IDs";
        for(i=0; i<numberSNPs; ++i)
            ofs << "#" << v.vecSNPNames[i];
        ofs << endl;

        unsigned long long segIndex(0);       // each char for 4 segments [ 11 01 00 10 ]
        const unsigned long long numx(numberSNPs*numberSamples);
        unsigned char shift;
        for(i=0; i<numberSamples; ++i)
        {
            if(i%7==0)  // show % every 7 samples
                progShow('W', unsigned(100*i/numberSamples));

            ofs << left << setw(32) << v.vecSampleIDs[i];
            for(j=0; j<numx; j+=numberSamples)
            {
                //              segIndex =  j * numberSamples   +   i;
                //          [base]   [offset]
                segIndex =  j    +   i;
                //                charIndex = segIndex / 4;           // the index of vecGenotype for the current Animal
                /*
                 // determine the position of genotype for the current animal
                 //     2bit  3*2 0s in the left    curSeg
                 shift= 2 *   (3                  - segIndex % 4);

                 mask =  3 << shift;    //  0000 0011 shift left to get the mask

                 ctmp = vecGenotypes[charIndex];
                 ctmp &= mask;
                 // shift ctmp right back to get the coded genotype
                 ctmp >>= shift;

                 // at last,
                 //          mask          char = 4-seg-genotypes        shift
                 //ctmp = ( (3 << shift) & vecGenotypes[charIndex] ) >> (2 * (3 - segIndex % 4));
                 ofs << vecUC2STR[ctmp] << " ";

                 *****************
                 *
                 *   for speed optimization
                 *   detailed algorithm see above
                 */
                shift = 2 * (3 - (segIndex % 4));
                ofs  << " " << vecUC2STR[( (3 << shift) & v.vecGenotypes[segIndex/4] ) >> shift];
            }
            ofs << endl;
        }
        progClear();

        ofs.close();
        v.clear();
    }
}

int main(const int argc, const char * argv[])
{
    libcbk::SNPipeline_VCF2AB(argc, argv);
    return 0;
}
