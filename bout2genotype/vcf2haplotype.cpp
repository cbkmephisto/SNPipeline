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
//  vcf2haplotype.cpp
//  libHailins
//
//  Created by Hailin SU on 3/13/14.
//
//  vcf2haplotype

#include "../appInfo.h"
#include "VCF_Record.h"
#include <vector>
#include <sstream>

using namespace std;
namespace libcbk
{
    class vcf2haplotype : public appInfo
    {
    public:
        vcf2haplotype(const int argc, const char * argv[]);
    private:
        string outputPrefix;
        void theMain(const int argc, const char * argv[]);
        void send2write(const vector<string> & vec2write, const unsigned short chr);
    };

    vcf2haplotype::vcf2haplotype(int argc, const char * argv[])
    {
        progName    = "vcf2haplotype";
        version     = "2014-03-13 12:23 CDT";
        created     = "2014-03-13";
        updated     = "N/A";
        descrip     = "Read beagle4 phased output vcf format, converting to haplotype info.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2013-02-02 init
         */
    }

    void vcf2haplotype::theMain(int argc, const char * argv[])
    {
        makeSure(argc >= 2, argv[0], "imputed/phased.*.vcf");

        for(unsigned i=1; i<argc; ++i)
        {
            VCF_File vf(argv[i]);
            if(!vf.open())
            {
                vf.close();
                E("ss", "Unable to open vcf file:", argv[i]);
                continue;
            }

            // ################## load the file into vecAllWindows
            I("sss", 1, "loading beagle 4 output vcf file : [", argv[i], "]");
            CQueue cqtmp;
            vector<CQueue> vecAllWindows;
            while((cqtmp=vf.getNextQueue()).size())
            {
                progShow('L', unsigned(vf.curTellG/vf.FILE_SIZE_1p));
                vecAllWindows.push_back(cqtmp);
            }
            progClear();
            unsigned j, k, z, lastChr(0);

            outputPrefix = argv[i];
            outputPrefix = "haplotype.from." + outputPrefix;

            I("s", 1, "writing haplotypes out");
            progShow('W', 0);
            CQueue * cq;
            VCF_Record * cr;
            vector<string> vec2write;
            lastChr=vecAllWindows[0].curChr;
            stringstream ss;
            string strChr, strWin;
            ss << lastChr;
            ss >> strChr;
            ss.clear();
            vec2write.resize(2*vf.vecSampleIDs.size());
            for(z=0; z<vec2write.size(); z++)
                vec2write[z]=vf.vecSampleIDs[z/2];// + "\tchr" + strChr + "\t";

            unsigned p2p(0.02*vecAllWindows.size());
            if(p2p==0)
                p2p=1;
            for(j=0; j<vecAllWindows.size(); ++j)   // for each window
            {
                // vecAllWindows[j] is a CQueue = window = vector<VCF_Record>
                cq=&vecAllWindows[j];
                if(lastChr!=cq->curChr) // new chr, new block
                {
                    send2write(vec2write, lastChr);
                    lastChr=cq->curChr;
                    ss << lastChr;
                    ss >> strChr;
                    ss.clear();
                    for(z=0; z<vec2write.size(); z++)
                    {
                        vec2write[z]=vf.vecSampleIDs[z/2];// + "\tchr" + strChr + "\t";
//                        vec2write[z].shrink_to_fit();
                    }
                }

                ss << cq->curWindow;
                ss >> strWin;
                ss.clear();
                for(z=0; z<vec2write.size(); z++)
                {
                    vec2write[z].push_back(' ');
                    vec2write[z] += strChr;
                    vec2write[z].push_back('_');
                    vec2write[z] += strWin;
                    vec2write[z].push_back('_');
                }
                for(k=0; k<cq->size(); ++k) // for each record
                {
                    cr=&(*cq)[k];
                    for(z=0; z<vec2write.size(); z++)
                        vec2write[z].push_back('0'+(*cr).vec_b_GenotypeLine[z]);
                }
                if(j%p2p==0)
                    progShow('W', j*100/vecAllWindows.size());
            }
            send2write(vec2write, lastChr);
            progClear();

            vector<string>().swap(vec2write);
            vf.close();

            for(j=0; j<vecAllWindows.size(); ++j)
                vector<VCF_Record>().swap(vecAllWindows[i]);
            vector<CQueue>().swap(vecAllWindows);
        }
    }


    void vcf2haplotype::send2write(const vector<string> &vec2write, const unsigned short curChr)
    {
        static unsigned i;
        static stringstream sst;
        static string sChr;
        sst.clear();
        sst << curChr;
        sst >> sChr;

        string ofn (outputPrefix + "." + sChr);
        ofstream ofs (ofn.c_str());
        if(!ofs)
        {
            cout << "Do not have the write access in current folder. Canceling." << endl;
            exit (-150);
        }

        for(i=0; i<vec2write.size(); ++i)
            ofs << vec2write[i] << endl;
        ofs.close();
    }

}

int main(int argc, const char * argv[])
{
    libcbk::vcf2haplotype(argc, argv);
    return 0;
}