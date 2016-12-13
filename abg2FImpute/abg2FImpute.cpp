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
//  abg2FImpute.cpp
//  libHailins
//
//  Created by Hailin SU on 3/27/14.
//

#include "abg2FImpute.h"
//  abg2FImpute

using namespace std;
namespace libcbk
{
    class abg2FImpute : public appInfo
    {
    public:
        abg2FImpute(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
    };

    abg2FImpute::abg2FImpute(int argc, const char * argv[])
    {
        progName    = "abg2FImpute";
        version     = "2014-10-30 18:20 CDT";
        created     = "2014-03-27";
        updated     = "2014-04-23";
        descrip     = "Proceed ab-genotype file to be a FImpute format. A snp_info.txt file is needed.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-03-27 init
         2014-04-23 corrected SNP order mismatch problem
         */
    }

    void abg2FImpute::theMain(int argc, const char * argv[])
    {
        makeSure(argc >= 3, argv[0], "snp_info.txt genotype1 (*)");

        //############################################# read the snp_info
        I("ss", 1, "reading", argv[1]);
        ifstream infs (argv[1]);
        if(!infs)
        {
            E("ss", "File not found:", argv[1]);
            exit(-1);
        }
        splitablestring aLine;
        stringstream sst;
        // read the header to determine how many chips
        getline(infs, aLine);
        sst << aLine;
        //      SNP_ID  Chr     Pos
        sst >> aLine >> aLine >> aLine;
        unsigned noc(0); // number of chips
        while(sst>>aLine)
            ++noc;
        sst.clear();

        vector< map<string, unsigned> > vecMaps;    // vector of mapSNPName2Pos
        vecMaps.resize(noc);

        string SNP_ID;
        unsigned i, utmp, linesR(0);
        unsigned long fs (0.01*tellFileSize(argv[1]));
        if(fs<1)
            fs=10;
        while(getline(infs, aLine))
        {
            sst << aLine;
            //               Chr        Pos
            sst >> SNP_ID >> aLine >> aLine;
            
/*
            for(i=0; i<noc; ++i)
            {
                sst >> utmp;    // the Chipx column of current line (SNP)
                if(utmp>0)      // if not 0, store the index of SNPName
                {
                    vecMaps[i][SNP_ID]=utmp-1;
//                    if(i==0 && utmp<10)
//                        cout << utmp-1 << " " << SNP_ID << endl;
                }
            }
*/
            i=0;
            while(sst>>utmp)
            {
                if(utmp>0)      // if not 0, store the index of SNPName
                    vecMaps[i][SNP_ID]=utmp-1;
                ++i;
            }
            sst.clear();

            if(++linesR%13000==0)
                progShow('L', unsigned(infs.tellg()/fs));
            
            if(i!=noc)
            {
                E("si", "bad map file: column number not match for line", 1+linesR);    // 1 for header line
                infs.close();
                exit (-2);
            }
                
        }
        infs.close();
        progClear();
        I("is", 2, noc, "chips read.");

        //############################################# process genotype files in order
        // test if the input files are readable, and reorder the input file by chip order
        unsigned long nos;   // number of snps
        vector<string> vecFn;
        vecFn.resize(noc);
        for(i=0; i<noc; ++i)
            vecFn[i]="";

        for(unsigned argvi=2; argvi<argc; ++argvi)  // for each genotype file
        {
            infs.open(argv[argvi]);
            if(!infs)
            {
                E("ss", "File not found:", argv[argvi]);
                exit (-2);
            }
            getline(infs, aLine);                   // read the first line to get the SNP number
            bool matches(false);
            nos=aLine.split('#').size()-1;
            for(i=0; i<noc; ++i)
                if(vecMaps[i].size()==nos)
                {
                    matches=true;
                    vecFn[i]=argv[argvi];
                }
            infs.close();
            if(!matches)
            {
                E("sisss", "The number of SNPs", nos, "in the given file", argv[argvi], "does not match any chip map in", argv[1]);
                for(i=0; i<noc; ++i)
                    I("sisiss", 2, "Chip", i+1, ".size =", vecMaps[i].size(), "--", vecFn[i].c_str());
                exit (-3);
            }
        }

        for(i=0; i<noc; ++i)
            I("sisiss", 2, "Chip", i+1, ".size =", vecMaps[i].size(), "--", vecFn[i].c_str());

        //############################################# go go go
        ofstream ofs ("genotypes.txt");
        if(!ofs)
        {
            infs.close();
            E("s", "could not open file to write: genotypes.txt");
            exit(-11);
        }
        ofs << setfill(' ') << left
        << setw(32) << "ID" << right
        << setw(6) << "Chip"
        << "    " << "Call..." << endl;

        vector<string> vt;  // vec tmp
        vector<unsigned> vecInd;
        vector<char> vecReorderedLine;
        map<string, string> mapCoding;
        mapCoding["AA"]="0";
        mapCoding["AB"]="1";
        mapCoding["BB"]="2";
        mapCoding["--"]="5";
        I("s", 1, "converting...");
        for(unsigned cn=0; cn<noc; ++cn)    // foreach input file, from higher density to lower one
        {
            if(vecFn[cn].length()==0)   // skip empty entry
                continue;
            I("s", 2, vecFn[cn].c_str());
            infs.open(vecFn[cn].c_str());
            unsigned long fs (0.01*tellFileSize(vecFn[cn]));
            if(fs<1)
                fs=10;
            // reorder from the SNP Names
            getline(infs, aLine);
//            I("s", 3, "splitting aLine to vt");
            vt=aLine.split('#');
//            I("s", 3, "resizing vecInd");
            vecInd.resize(vt.size()-1);
//            I("s", 3, "resizing vecReorderedLine");
            vecReorderedLine.resize(vt.size()-1);
//            I("s", 3, "filling vecInd");
            for(i=1; i<vt.size(); ++i)  // calc the new index of each column
            {
//                vecInd[i-1]=vecMaps[cn][vt[i]];
                vecInd[vecMaps[cn][vt[i]]]=i-1;
//                if(vecMaps[cn][vt[i]]<8) // print out the first 8 abgs
//                {
//                    cout << "vecInd[" << vecMaps[cn][vt[i]] << "] = " << i-1 << " -> " << vt[i] << endl;
//                }
            }
            linesR=0;
            string *spt;
            
//            I("s", 3, "writing");
            while(getline(infs, aLine)) // foreach line
            {
                vt=aLine.split();
                if(vt[0].length()>29)   // skip if ID_length>29
                    continue;
                ofs << setfill(' ') << left
                << setw(32) << vt[0] << right   // ID
                << setw(6) << cn+1              // chip
                << "    ";
                aLine=string("");
                for(i=0; i<vt.size()-1; ++i)  // reordering
                {
//                    SOMETHING WRONG AROUND HERE!!!!!!!!!!!!!!!!!!!!!
                    spt=& (vt[1+vecInd[i]]);
//                    if(linesR<2 && i<9) // print out the first 8 abgs
//                    {
//                        cout << i << " -> vecInd[" << i << "] = " << vecInd[i] << " -> " << (*spt) << endl;
//                    }
                    aLine+=  (mapCoding.find(*spt)==mapCoding.end()? "5" : mapCoding[*spt]);    // get all missing
                }
                ofs << aLine << endl;
                if(++linesR%100==1)
                    progShow('W', unsigned(infs.tellg()/fs));
            }// end of foreach line
            infs.close();
            progClear();

//            I("s", 3, "freeing vt");
            vector<string>().swap(vt);
//            I("s", 3, "freeing vecInd");
            vector<unsigned>().swap(vecInd);
//            I("s", 3, "freeing vecReorderedLine");
            vector<char>().swap(vecReorderedLine);
        }// end of foreach input file
        ofs.close();

        // generating dummy ctr file
        ofs.open("dummy.ctr");
        ofs
        << "title=\"abg2FImpute dummy ctr\";" << endl
        << "genotype_file=\"genotypes.txt\";" << endl
        << "snp_info_file=\"" << argv[1] << "\" /chrx=30;" << endl
        << "exclude_chr=33;" << endl
        << "output_folder=\"optDummy\";" << endl
        << "njob=8;" << endl
        << "ped_file=\"pedFImpute\";" << endl
        << "parentage_test /remove_conflict;" << endl;
        
        ofs.close();
    }
}

int main(int argc, const char * argv[])
{
    libcbk::abg2FImpute(argc, argv);
    return 0;
}