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
//  fout2haplotype.cpp
//  libHailins
//
//  Created by Hailin SU on 4/21/14.
//

#include "fout2haplotype.h"

//  fout2haplotype

using namespace std;
namespace libcbk
{
    class fout2haplotype : public appInfo
    {
    public:
        fout2haplotype(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
    };

    fout2haplotype::fout2haplotype(int argc, const char * argv[])
    {
        progName    = "fout2haplotype";
        version     = "2016-04-25 16:14 CDT";
        created     = "2014-04-21";
        updated     = "2016-04-25";
        descrip     = "Converts FImpute output file to haplotype file.\n\t+/- in the output SNP loci means (heterozygous + unable to phase)\n\t X  in the output SNP loci means unable to impute. Will read a file windowSize for window size if exists. Default 1M.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-04-21 init
         2015-04-07 dealing with output still unimputed code '5'
         2015-06-04 filling the gap between two adjacent windows
         2016-04-25 windowSize
         */
    }

    void fout2haplotype::theMain(int argc, const char * argv[])
    {
        makeSure(argc == 5 || argc == 7 , argv[0], "snp_info.txt genotypes_imp.txt -o haplotype.from.fout [ chr win ]");
        makeSure(argv[3] == string("-o"), argv[0], "snp_info.txt genotypes_imp.txt -o haplotype.from.fout [ chr win ]");

        string chr, win;
        if(argc==7)
        {
            chr=argv[5];
            win=argv[6];
        }
        I("s", 1, "Reading snp_info.txt ...");

        VecWindows vw(argv[1]);

        ifstream infs(argv[2]);
        if(!infs)
        {
            E("ss", "File not found:", argv[2]);
            exit(-2);
        }

        ofstream ofs(argv[4]);
        if(!ofs)
        {
            infs.close();
            E("ss", "Could not open file to write:", argv[4]);
            exit(-3);
        }

        unsigned long long f1p(0.01*tellFileSize(argv[2]));
        if(f1p<1)
            f1p=1;

        unsigned linesR(0);
        string stmp;

        getline(infs, stmp);    // get rid of header line

        I("s", 1, "Converting genotypes_imp.txt to haplotype file...");
        while(getline(infs, stmp))
        {
            if(argc==7)
                vw.load(stmp, chr, win);
            else
                vw.load(stmp);
            ofs << vw.lineP << endl;
            ofs << vw.lineM << endl;
            if(++linesR%200==1)
                progShow('P', infs.tellg()/f1p);
        }
        progClear();

        I("ss", 1, argv[4], "generated.");
        infs.close();
        ofs.close();
    }
}

int main(int argc, const char * argv[])
{
    libcbk::fout2haplotype(argc, argv);
    return 0;
}