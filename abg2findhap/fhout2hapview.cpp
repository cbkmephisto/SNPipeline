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
//  fhout2plink.cpp
//  libHailins
//
//  Created by Hailin Su on 2015-08-10.
//

#include "../appInfo.h"
#include "../splitablestring.h"
#include <map>
#include <iomanip>
#include <set>

using namespace std;
namespace libcbk
{
    class fhout2hapview : public appInfo
    {
    public:
        fhout2hapview(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
        bool parseParameters(const int argc, const char* argv[]);
        string fnP, fnM, fnG, fnO, fnI;
        unsigned tgtChr, sttWin, nWin;
    };

    fhout2hapview::fhout2hapview(int argc, const char * argv[])
    {
        progName    = "fhout2hapview";
        version     = "2015-08-10 16:15 CDT";
        created     = "2015-08-10";
        updated     = "N/A";
        descrip     = "Converts findhap output to HapView input to plot LD...";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2015-08-10 init
         2015-08-14 add -i, -w, -n options
         */
    }
    
    void fhout2hapview::theMain(int argc, const char * argv[])
    {
        makeSure(parseParameters(argc, argv), argv[0], "[-i inc] -c 14 -w 24 -n 1 -p pedigree.file -m chromosome.data -g haplotypes.txt -o PREFIX");
        
        //############################# convert chromosome.dat to PREFIX.map
        splitablestring aLine;
        vector<string> vtmp;

        ifstream infs (fnM.c_str());
        if(!infs)
        {
            E("ss", "input file fnM not found:", fnM.c_str());
            exit(-20);
        }
        
        string fout (fnO + ".map");
        ofstream ofs (fout.c_str());
        if(!ofs)
        {
            infs.close();
            E("ss", "could not open file to write:", fout.c_str());
            exit(-21);
        }
        I("ssss", 1, "converting", fnM.c_str(), "to", fout.c_str());
        // skip header
        getline(infs, aLine);
        // keep 2[chr] 1[rs] '0'[genetic distance] 5[pos]
        unsigned i, curChr;
        int begInd(-1), endInd(-1), curInd(-1);
        string snpId;
        unsigned long curPos;
        stringstream sst;
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            if(vtmp.size()<5)
            {
                E("sssiss", "not a 5+col through file:", fnM.c_str(), "-", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                ofs.close();
                exit(-20);
            }
            sst << aLine;
            //     snpName  chr       within    overall   location       ;  ignore the rest
            sst >> snpId >> curChr >> curPos >> curPos >> curPos;
            sst.str("");
            sst.clear();

            ++curInd;

            // test chr
            if(tgtChr!=curChr && begInd<0)  // before the tgtChr
                continue;
            
            // test pos, 20150814
            if(sttWin*1000000 > curPos && begInd<0)  // before the starting win
                continue;

            if((sttWin+nWin)*1000000 < curPos && begInd)  // after n win
                break;

            if(begInd<0)
                begInd=curInd;

            ofs << setfill(' ') << right
//            << setw(8)  << vtmp[1]      // chr
            << setw(32) << vtmp[0]      // ID
//            << setw(8)  << 0            // genetic distance
            << setw(12) << vtmp[4]      // pos
            << endl;
        }
        endInd=curInd-1;            // inclusive endInd!
        infs.close();
        ofs.close();
        
        //############################# convert pedigree.file and haplotypes.txt to PREFIX.ped
        // first read include file
        I("s", 1, "converting ped file");
        I("s", 2, "deal with inc file (-i)");
        set<string> setInc;
        if(fnI!="")
        {
            infs.open(fnI.c_str());
            if(!infs)
                E("sss", "input file not found:", fnI.c_str(), " IGNORING INC OPTION (-i)");
            else
            {
                while(getline(infs, aLine))
                    setInc.insert(aLine);
            }
            infs.close();
        }
        // read pedigree.file
        I("ss", 2, "loading", fnP.c_str());
        map<string, string> mapID2Gender;
        map<string, string> mapID2Sire;
        map<string, string> mapID2Dam;
        infs.open(fnP.c_str());
        if(!infs)
        {
            E("ss", "input file not found:", fnP.c_str());
            exit(-20);
        }
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            if(vtmp.size()!=7)
            {
                E("sssiss", "not a 7 col through file:", fnP.c_str(), " - ", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                ofs.close();
                exit(-20);
            }
            // vtmp[1] = ID
            //      2    FID
            //      3    MID
            //      0    SEX
            mapID2Gender[vtmp[1]]=vtmp[0];
            mapID2Sire[vtmp[1]]  =vtmp[2];
            mapID2Dam[vtmp[1]]   =vtmp[3];
        }
        infs.close();
        
        // and then read haplotype.txt and write out the result
        infs.open(fnG.c_str());
        if(!infs)
        {
            E("ss", "input file not found:", fnG.c_str());
            exit(-20);
        }
        
        fout=fnO+".ped";
        ofs.open(fout.c_str());
        if(!ofs)
        {
            infs.close();
            E("ss", "could not open file to write:", fout.c_str());
            exit(-21);
        }
        I("ssss", 2, "converting", fnG.c_str(), "to", fout.c_str());
        map<string, unsigned> mapGender;
        mapGender["M"]=1;
        mapGender["F"]=2;
        mapGender["U"]=0;
        unsigned vtmp3length;
        set<string> setGenotypedIDs;
        string famID, sireID, damID;
        unsigned uFamID(0);
        while(getline(infs, aLine))
        {
            vtmp=aLine.split();
            // [0] ID
            // [1] chip
            // [2] length
            // [3] genotype
            if(vtmp.size()!=4)
            {
                E("sssiss", "not a 4-col through file:", fnG.c_str(), "-", vtmp.size(), "fields found in line starting with", vtmp[0].c_str());
                infs.close();
                ofs.close();
                exit(-20);
            }
            sireID=mapID2Sire[vtmp[0]];
            if(sireID[0]=='-')
                sireID="0";
            damID=mapID2Dam[vtmp[0]];
            if(damID[0]=='-')
                damID="0";
            if(sireID==damID && sireID=="0")
            {
                sst << ++uFamID;
                sst >> famID;
                sst.clear();
            }
            else
                famID=mapID2Sire[vtmp[0]]+"-"+ mapID2Dam[vtmp[0]];
            
            if(!setInc.size() || setInc.find(vtmp[0])!=setInc.end())    // only output if setInc is empty, or ID in setInc
            {
                ofs << setfill(' ') << left
                << setw(32) << famID << right     // family ID
                << setw(16) << vtmp[0]                                                  // ID
                << setw(16) << sireID                                      // sireID
                << setw(16) << damID                                       // damID
                << setw(16) << mapGender[mapID2Gender[vtmp[0]]]                         // sex
                << " 1 ";                                                              // phe
                aLine=string("");
                vtmp3length=(unsigned)vtmp[3].length();
                for(i=begInd*2; i<=endInd*2+1 && i<vtmp3length; ++i)
                {
                    aLine+=" ";
                    aLine+=vtmp[3][i];
                    if(i%2)
                        aLine+=" ";
                }
                ofs << aLine << endl;
                setGenotypedIDs.insert(vtmp[0]);
            }
        }
        infs.close();

        // not-genotyped: set to missing phe
        map<string, string>::iterator mit;
        for(mit=mapID2Sire.begin(); mit!=mapID2Sire.end(); ++mit)
        {
            aLine=mit->first;
            if(setGenotypedIDs.find(aLine)==setGenotypedIDs.end()) // not genotyped ancestor
            {
                sireID=mapID2Sire[aLine];
                if(sireID[0]=='-')
                    sireID="0";
                damID=mapID2Dam[aLine];
                if(damID[0]=='-')
                    damID="0";
                if(sireID==damID && sireID=="0")
                {
                    sst << ++uFamID;
                    sst >> famID;
                    sst.clear();
                }
                else
                    famID=mapID2Sire[aLine]+"-"+ mapID2Dam[aLine];

                if(!setInc.size() || setInc.find(vtmp[0])!=setInc.end())    // only output if setInc is empty, or ID in setInc
                {
                    
                    ofs << setfill(' ') << left
                    << setw(32) << famID << right     // family ID
                    << setw(16) << aLine                                                  // ID
                    << setw(16) << sireID                                      // sireID
                    << setw(16) << damID                                       // damID
                    << setw(16) << mapGender[mapID2Gender[aLine]]                         // sex
                    << " -9 ";                                                              // phe, missing
                    
                    aLine=string("");
                    for(i=begInd*2; i<=endInd*2+1; ++i)                         // genotype, missing
                    {
                        aLine+=" 0";
                        if(i%2)
                            aLine+=" ";
                    }
                    ofs << aLine << endl;
                }
            }
        }
        ofs.close();
    }
    
    bool fhout2hapview::parseParameters(const int argc, const char* argv[])
    {
        stringstream sst;
        string stmp;
        bool ret(true);
        fnP=fnM=fnG=fnO=fnI="";
        tgtChr=sttWin=nWin=0;
        for(unsigned i=1; i<argc; ++i)
        {
            stmp=argv[i];
            if(stmp=="-c") // read the target Chr
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-c" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                {
                    stmp=argv[++i];     // convert it to unsigned
                    sst << stmp;
                    sst >> tgtChr;
                    sst.str("");
                    sst.clear();
                }
            }
            else if(stmp=="-w") // starting win
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-w" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                {
                    stmp=argv[++i];     // convert it to unsigned
                    sst << stmp;
                    sst >> sttWin;
                    sst.str("");
                    sst.clear();
                }
            }
            else if(stmp=="-n") // number windows
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-n" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                {
                    stmp=argv[++i];     // convert it to unsigned
                    sst << stmp;
                    sst >> nWin;
                    sst.str("");
                    sst.clear();
                }
            }
            else if(stmp=="-i") // read the inc file
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ;   // do nothing; just ignore
                else
                    fnI=argv[++i];     // stores the include-fileName and open it for reading
            }
            else if(stmp=="-p") // read the pedigree files
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
        if(!tgtChr)
        {
            E("s", "ERROR: target chr (-c) should be specified.");
            ret=false;
        }
        if(!sttWin)
        {
            E("s", "ERROR: starting win (-w) should be specified.");
            ret=false;
        }
        if(!nWin)
        {
            E("s", "ERROR: number of windows (-n) should be specified.");
            ret=false;
        }
        return ret;
    }
}

int main(int argc, const char * argv[])
{
    libcbk::fhout2hapview(argc, argv);
    return 0;
}