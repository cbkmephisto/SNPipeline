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
//  fhout2fiout.cpp
//  libHailins
//
//  Created by Hailin SU on 2014-11-16.
//

//===== replacing area =====
//  cpp filename: fhout2fiout
//  date        : 2014-11-16
//  time        : 11:55

#include "../appInfo.h"
#include "../splitablestring.h"
#include <map>
#include <iomanip>

using namespace std;
namespace libcbk
{
    class fhout2fiout : public appInfo
    {
    public:
        fhout2fiout(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
        bool parseParameters(const int argc, const char* argv[]);
        string fnP, fnM, fnG, fnO;
    };
    
    fhout2fiout::fhout2fiout(int argc, const char * argv[])
    {
        progName    = "fhout2fiout";
        version     = "2014-11-16 11:55 CDT";
        created     = "2014-11-16";
        updated     = "2014-11-16";
        descrip     = "extracts things like this";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-11-16 init
         */
    }
    
    void fhout2fiout::theMain(int argc, const char * argv[])
    {
        makeSure(parseParameters(argc, argv), argv[0], "-p pedigree.file -m chromosome.data -g haplotypes.txt -o PREFIX");
        
        //############################# read animalID xref
        I("sss", 1, "reading", fnP.c_str(), "for animalID xref");
        ifstream infs (fnP.c_str());
        if(!infs)
        {
            E("ss", "input file not found:", fnP.c_str());
            exit(-20);
        }
   
        // [1] -> [6]
        map<string, string> mapiID2IID;
        splitablestring aLine;
        vector<string> vtmp;
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            if(vtmp.size()!=7)
            {
                E("sssiss", "not a 7-col through file:", fnP.c_str(), "-", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                exit(-20);
            }
            mapiID2IID[vtmp[1]]=vtmp[6];
        }
        infs.close();
        vector<string>().swap(vtmp);

        //############################# convert chromosome.dat to snp_info.txt
        infs.open(fnM.c_str());
        if(!infs)
        {
            E("ss", "input file not found:", fnM.c_str());
            exit(-20);
        }

        string fout (fnO + ".snp_info.txt");
        ofstream ofs (fout.c_str());
        if(!ofs)
        {
            infs.close();
            E("ss", "could not open file to write:", fout.c_str());
            exit(-21);
        }
        I("ssss", 1, "converting", fnM.c_str(), "to", fout.c_str());
        // keep 0 1 4 6 +
        unsigned i;
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            if(vtmp.size()<7)
            {
                E("sssiss", "not a 7+col through file:", fnM.c_str(), "-", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                ofs.close();
                exit(-20);
            }
            // TODO: copy mapUniter how to make snp_info.txt
            ofs << setfill(' ') << left
            << setw(55) << vtmp[0] << right
            << setw(5)  << vtmp[1]
            << setw(12) << vtmp[4];
            for(i=6; i<vtmp.size(); ++i)
                ofs << setfill(' ') << right << setw(14) << vtmp[i];
            ofs << endl;

        }
        infs.close();
        ofs.close();
        
        //############################# convert haplotypes.txt to genotypes_imp.txt
        infs.open(fnG.c_str());
        if(!infs)
        {
            E("ss", "input file not found:", fnG.c_str());
            exit(-20);
        }
        
        fout=fnO+".genotypes_imp.txt";
        ofs.open(fout.c_str());
        if(!ofs)
        {
            infs.close();
            E("ss", "could not open file to write:", fout.c_str());
            exit(-21);
        }
        I("ssss", 1, "converting", fnG.c_str(), "to", fout.c_str());
        ofs << setfill(' ') << left
        << setw(30) << "ID" << right
        << setw(5) << "Chip"
        << " " << "Call..." << endl;
        // keep xref([0]) 1 3
        map<unsigned short, string> mapGenCoder;
        // I know it is hard to understand
        // 1A 2B
        // "11" = '1' '1' = 49 49 = 0x31 0x31 -->> 0x3131 -> AA=0
        // "12" = '1' '2' = 49 50 = 0x31 0x32 -->> 0x3231 -> AB=3   0x31 0x32 => 0x3231 is determined by the way the computer stores data
        // "21" = '2' '1' = 50 49 = 0x32 0x31 -->> 0x3132 -> BA=4   it is not a mistake here
        // "22" = '2' '2' = 50 50 = 0x32 0x32 -->> 0x3232 -> BB=2
        mapGenCoder[0x3131]="0";
        mapGenCoder[0x3132]="4";
        mapGenCoder[0x3231]="3";
        mapGenCoder[0x3232]="2";
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            if(vtmp.size()!=4)
            {
                E("sssiss", "not a 4-col through file:", fnG.c_str(), "-", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                ofs.close();
                exit(-20);
            }
            ofs << setfill(' ') << left
            << setw(30) << mapiID2IID[vtmp[0]] << right   // ID
            << setw(5)  << vtmp[1]              // chip
            << " ";
            aLine=string("");
            for(i=0; i<vtmp[3].length(); i+=2)
                // I know it is hart to understand
                aLine+=mapGenCoder[*(short*)&(vtmp[3][i])];
            ofs << aLine << endl;
        }
        infs.close();
        ofs.close();
    }
    
    bool fhout2fiout::parseParameters(const int argc, const char* argv[])
    {
        string stmp;
        bool ret(true);
        fnP=fnM=fnG=fnO="";
        for(unsigned i=1; i<argc; ++i)
        {
            stmp=argv[i];
            if(stmp=="-p") // read the pedigree files
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-p" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                    fnP=argv[++i];     // stores the list-fileName and open it for reading
            }
            else if(stmp=="-m") // inserts the mapFile
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                    fnM=argv[++i];
            }
            else if(stmp=="-g") // inserts the genFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                    fnG=argv[++i];
            }
            else if(stmp=="-o") // inserts the genFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                    fnO=argv[++i];
            }
            else
                ret=false;
        }
        if(fnP.length()==0)
        {
            E("s", "ERROR: No ped file (-p) specified.");
            ret=false;
        }
        if(fnM.length()==0)
        {
            E("s", "ERROR: No map file (-m) specified.");
            ret=false;
        }
        if(fnG.length()==0)
        {
            E("s", "ERROR: No gen file (-g) specified.");
            ret=false;
        }
        if(fnO.length()==0)
        {
            E("s", "ERROR: No out file (-o) specified.");
            ret=false;
        }
        return ret;
    }
}

int main(int argc, const char * argv[])
{
    libcbk::fhout2fiout(argc, argv);
    return 0;
}