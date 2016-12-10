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
//  abg2findhap.cpp
//  libHailins
//
//  Created by Hailin SU on 2014-11-10.
//  Copyright (c) 2013 iastate. All rights reserved.
//

//===== replacing area =====
//  cpp filename: abg2findhap
//  date        : 2014-11-10

#include "../appInfo.h"
#include "abg2findhap.h"

using namespace std;
namespace libcbk
{
    //*********************************************************************************************
    //************************     abg2findhap     ************************************************
    //*********************************************************************************************
    abg2findhap::abg2findhap(int argc, const char * argv[])
    {
        progName    = "abg2findhap";
        version     = "2014-11-10 10:04 CDT";
        created     = "2014-11-10";
        updated     = "2014-11-10";
        descrip     = "Convert ab-genotype files with pedigree and map info to be imputable by findhap V3.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2014-11-10 init
         */
    }
    
    void abg2findhap::readme()
    {
        cout << "\n\
        ======================================================================\n\
        ========================        ReadMe        ========================\n\
        ======================================================================\n\
        \n\
        * usage\n\
        abg2findhap [-s] -p PEDIGREE -m MAP_FILE(s) -g GEN_FILE(s)\n\
        \n\
        * THE OPTIONAL -s OPTION\n\
        - generate output files by chromosome\n\
        - is off by default\n\
        \n\
        * THE PEDIGREE\n\
        - reads pedRefiner output .csv format(coma-separated-values) 3-col pedigree file, convert to findhap V3 format\n\
        and make xref for IDs\n\
        - no checking ability was implemented, so\n\
        - pedRefiner version >= 2014-11-10 is required so that the output is checked and sorted\n\
        - date of birth will be made up\n\
        - unknown gender will be converted to M\n\
        \n\
        * THE MAP(s)\n\
        - reads 'SNP_NAME CHR POS' 3-col-through, space/tab delimetered file\n\
        - skip 1st line as header\n\
        \n\
        * THE GEN_FILE(s)\n\
        - ab-genotype file from SNPipeline" << endl;
        exit(-1);
    }
    
    void abg2findhap::theMain(const int argc, const char * argv[])
    {
        bOutputByChr=false;
        
        if(!parseParameters(argc, argv))
            readme();
        
        PedConverter p(pedFileName);
        
        MapConverter m(vecMapFileName, bOutputByChr);
        
        GenConverter g(vecGenFileName, p, m, bOutputByChr);
    }
    
    bool abg2findhap::parseParameters(const int argc, const char* argv[])
    {
        string stmp;
        bool ret(true);
        pedFileName="";
        for(unsigned i=1; i<argc; ++i)
        {
            stmp=argv[i];
            if(stmp=="-s")      // seperate output by chr is on
            {
                bOutputByChr=true;
            }
            else if(stmp=="-p") // read the pedigree files
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-l" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                    pedFileName=argv[++i];     // stores the list-fileName and open it for reading
            }
            else if(stmp=="-m") // inserts the mapFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                {
                    for(++i; i<argc && argv[i][0]!='-'; ++i)
                        vecMapFileName.push_back(argv[i]);
                    if(i!=argc)         // if not the end of param list, put the index back 1 position.
                        --i;
                }
            }
            else if(stmp=="-g") // inserts the genFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                    ret=false;
                else
                {
                    for(++i; i<argc && argv[i][0]!='-'; ++i)
                        vecGenFileName.push_back(argv[i]);
                    if(i!=argc)         // if not the end of param list, put the index back 1 position.
                        --i;
                }
            }
            else
                ret=false;
        }
        if(pedFileName.length()==0)
        {
            E("s", "ERROR: No ped file    (-p) specified.");
            ret=false;
        }
        if(vecMapFileName.size()==0)
        {
            E("s", "ERROR: No map file(s) (-m) specified.");
            ret=false;
        }
        if(vecGenFileName.size()==0)
        {
            E("s", "ERROR: No gen file(s) (-g) specified.");
            ret=false;
        }
        return ret;
    }
    
    //*********************************************************************************************
    //************************     PedConverter     ***********************************************
    //*********************************************************************************************
    PedConverter::PedConverter(const string fileName)
    {
        //////// load
        I("s", 1, "Entering PedConverter");
        ifstream infs (fileName.c_str());
        if(!infs)
        {
            E("sss", "error trying to open", fileName.c_str(), "to read, aborting.");
            exit(-1);
        }
        ofstream ofs ("pedigree.file");
        if(!ofs)
        {
            infs.close();
            E("s", "error trying to open 'pedigree.file' to write, aborting.");
            exit(-1);
        }
        I("s", 2, "loading");
        map<string, char>       mapiID2Gender;
        vector<string>          veciID;
        vector<string>          vecsID;
        vector<string>          vecdID;
        splitablestring         aLine;
        vector<string>          vecTmp;
        bool                    valid(true);
        unsigned                i(1);
        while(getline(infs, aLine))
        {
            vecTmp=aLine.split(',',1);
            if(vecTmp.size()==3)
            {
                veciID.push_back(vecTmp[0]);
                vecsID.push_back(vecTmp[1]);
                vecdID.push_back(vecTmp[2]);
                mapiID2nID[vecTmp[0]]=i++;
                // check predifened gender by pedigree: sire/dam
                if(vecTmp[1]!="0")
                {
                    if(mapiID2Gender.find(vecTmp[1])==mapiID2Gender.end())
                        mapiID2Gender[vecTmp[1]]='M';
                    else if(mapiID2Gender[vecTmp[1]]!='M')  // conflict
                    {
                        E("sssiss", "ERROR -", vecTmp[1].c_str(), "was NOT a MALE but listed as sire here, line", i, "starts with", vecTmp[0].c_str());
                        valid=false;
                    }
                }
                
                if(vecTmp[2]!="0")
                {
                    if(mapiID2Gender.find(vecTmp[2])==mapiID2Gender.end())
                        mapiID2Gender[vecTmp[2]]='F';
                    else if(mapiID2Gender[vecTmp[2]]!='F')  // conflict
                    {
                        E("sssiss", "ERROR -", vecTmp[2].c_str(), "was NOT a FEMALE but listed as dam here, line", i, "starts with", vecTmp[0].c_str());
                        valid=false;
                    }
                }
            }
            else
            {
                valid=false;
                E("siss", "Error reading ped file: not a 3-col csv file, line", i+1, "starting with", vecTmp[0].c_str());
                E("si", "    number of columns read:", vecTmp.size());
            }
        }
        infs.close();
        
        if(!valid)
        {
            exit(-2);
        }
        vector<string>().swap(vecTmp);
        
        //////// write
        I("s", 2, "writing");
        char gender;
        unsigned short bdy(2014);
        unsigned short bdm(11);
        unsigned short bdd(10);
        set<string> setUnknown;
        string aid;
        int nsID, ndID;
        for(i=0; i<veciID.size(); ++i)
        {
            aid=veciID[i];
            // fill gender
            if(mapiID2Gender.find(aid)==mapiID2Gender.end())  // info not found
            {
                if(aid.length()<6)
                {
                    setUnknown.insert(aid);
                    gender='M';
                }
                else
                    gender=aid[6];
                if(gender!='M' && gender!='F')
                {
                    setUnknown.insert(aid);
                    gender='M';
                }
            }
            else
                gender=mapiID2Gender[aid];
            
            // write the line
            nsID=mapiID2nID[vecsID[i]]>0? mapiID2nID[vecsID[i]]: -1-(i%2);
            ndID=mapiID2nID[vecdID[i]]>0? mapiID2nID[vecdID[i]]: -2+(i%2);
            ofs << setfill(' ') << left
            << setw(4)  << gender
            << right
            << setw(9) << mapiID2nID[aid]
            << " "
            << setw(9) << nsID
            << " "
            << setw(9) << ndID
            << " "
            << setw(5)  << bdy
            << setfill('0')
            << setw(2)  << bdm
            << setw(2)  << bdd
            << " "          // in case of too-long an animal ID
            << setfill(' ')// << "        " << left
            << setw(17) << mapiID2nID[aid]
            << " "
            << setw(30) << aid << endl;
            
            // recalc date
            if(++bdd>28)
            {
                if(++bdm>12)
                {
                    ++bdy;
                    bdm=1;
                }
                bdd=1;
            }
        }
        ofs.close();
        
        if(setUnknown.size()!=0)
        {
            E("s", "warning: unable to detect gender for (assigned to M arbitrarily)");
            i=0;
            for(set<string>::iterator sit=setUnknown.begin(); sit!=setUnknown.end(); ++sit)
            {
                cerr << "\t" << (*sit);
                if(++i%5==0)
                    cerr << endl;
            }
            if(i%5!=0)
                cerr << endl;
            set<string>().swap(setUnknown);
        }
        //////// clean
        map<string, char>().swap(mapiID2Gender);
        vector<string>().swap(veciID);
        vector<string>().swap(vecsID);
        vector<string>().swap(vecdID);
        
    }
    //*********************************************************************************************
    //************************     MapConverter     ***********************************************
    //*********************************************************************************************
    bool SNPMap::load(const char *fileName)
    {
        ifstream infs(fileName);
        if(!infs)
            return false;
        stringstream ss;
        splitablestring aLine;
        unsigned chr;
        unsigned pos;
        char idStr[12];   // cc ppp ppp ppp
        vector<string> vectmp;
        getline(infs, aLine); // skip the header line
        
        while(getline(infs, aLine))
        {
            vectmp=aLine.split();
            if(vectmp.size()!=3)
            {
                infs.close();
                cerr << "        Error reading map file: not a 3-col file" << endl
                << "            - line starting with:     " << vectmp[0] << endl
                << "            - number of columns read: " << vectmp.size() << endl;
                return false;
            }
            if(vectmp[1]=="X")
                vectmp[1]="30";
            else if(vectmp[1]=="MT")
                vectmp[1]="31";
            else if(vectmp[1]=="Y")
                vectmp[1]="33";
            ss << vectmp[1];
            ss >> chr;
            ss.str("");
            ss.clear();
            
            ss << vectmp[2];
            ss >> pos;
            ss.str("");
            ss.clear();
            
            if(chr==0 || pos==0)        // skip bad info
                continue;
            
            if(chr>30)                  // ignore chr>30!
                continue;
            
            sprintf(idStr, "%02u%09u", chr, pos);
            idStr[11]=0;
            aLine=idStr;
            setSNPIDs.insert(aLine);
            
            if(mapName2ID.find(vectmp[0])==mapName2ID.end())
                mapName2ID[vectmp[0]]=aLine;
            else if(mapName2ID[vectmp[0]]==aLine)
                cerr << "        WARNING(skipped): duplicated but same entry found for SNP " << vectmp[0] << endl;
            else
                mapName2ID[vectmp[0]]=aLine;    // multiple entries allowed in nameSNP->pos
            //                cerr << "        WARNING(skipped): duplicated SNPname entry with DIFFERENT position value found for SNP " << vectmp[0] << endl;
        }
        
        infs.close();
        return true;
    }
    
    bool SNPMap::containsID(const string &s)
    {
        return (setSNPIDs.find(s)!=setSNPIDs.end());
    }
    
    bool SNPMap::containsName(const string &s)
    {
        return (mapName2ID.find(s)!=mapName2ID.end());
    }
    
    unsigned SNPMap::sizeSetIDs()
    {
        return (unsigned)setSNPIDs.size();
    }
    
    unsigned SNPMap::sizeMapNames()
    {
        return (unsigned)mapName2ID.size();
    }
    
    // for highest density chip, adding alias in mapName2ID
    void SNPMap::mergeEntry(SNPMap &s)
    {
        static map<string, string>::iterator mit;
        for(mit=s.mapName2ID.begin(); mit!=s.mapName2ID.end(); ++mit) // for each name entry in ld chip
        {
            if(setSNPIDs.find(mit->second)!=setSNPIDs.end() && mapName2ID.find(mit->first)==mapName2ID.end()) // if ID in common and name is new
                mapName2ID[mit->first]=mit->second;
        }
    }
    
    // for lower density chips, deleting not-in-common SNPs
    void SNPMap::trimTo(SNPMap &s)
    {
        set<string> *psetValidSNPID (&(s.setSNPIDs));
        map<string, string> *pmapValidSNPName (&(s.mapName2ID));
        static set<string>::iterator sit;
        static map<string, string>::iterator mit;
        set<string> setNewSNPIDs;
        for(sit=setSNPIDs.begin(); sit!=setSNPIDs.end(); ++sit)   // for each SNP_ID in the current ld map
            if(psetValidSNPID->find(*sit)!=psetValidSNPID->end())   // if current SNP_ID in HD map, add to newSet
                setNewSNPIDs.insert(*sit);
        
        for(mit=mapName2ID.begin(); mit!=mapName2ID.end(); ++mit)   // for each SNP_Name in the current ld map
            if(pmapValidSNPName->find(mit->first)!=pmapValidSNPName->end())   // if current SNP_Name in HD map, add it's (HD version)ID to newSet
            {
                setNewSNPIDs.insert((*pmapValidSNPName)[mit->first]);
                if((*pmapValidSNPName)[mit->first]!=mit->second)
                {
                    cerr << "        WARNING(HD kept): " << mit->first << " is HD[" << (*pmapValidSNPName)[mit->first]
                    << "] and LD[" << mit->second << "]" << endl;
                }
            }
        setSNPIDs.swap(setNewSNPIDs);
        set<string>().swap(setNewSNPIDs);
    }
    
    MapConverter::MapConverter(vector<string> &vecFiles)
    {
        MapConverter(vecFiles, false);
    }
    
    MapConverter::MapConverter(vector<string> &vecFiles, bool bOutputByChr)
    {
        vector<SNPMap> vecSNPMapsX;
        I("s", 1, "Entering MapConverter");
        //////// load
        I("s", 2, "loading");
        map<unsigned, unsigned> mapN2iSNPMap;
        unsigned nMax(0);
        int ntmp;
        unsigned i;
        for(i=0; i<vecFiles.size(); ++i)
        {
            cout << "                    " << left << setfill(' ') << setw(16) << vecFiles[i].c_str() << flush;
            vecSNPMapsX.push_back(SNPMap());
            if(!vecSNPMapsX[i].load(vecFiles[i].c_str()))
            {
                cout << endl;
                E("sss", "error trying to open", vecFiles[i].c_str(), "to read, aborting.");
                exit(-2);
            }
            ntmp=vecSNPMapsX[i].sizeSetIDs();
            cout << right << setw(14) << ntmp << " markers read" << endl;
            mapN2iSNPMap[ntmp]=i;
            if(ntmp>nMax)
                nMax=ntmp;
        }
        vecChips.push_back(vecSNPMapsX[mapN2iSNPMap[nMax]]);
        
        //////// merge
        I("s", 2, "merging maps: adding nameSNP alias");
        i=2;
        unsigned delta(0);
        for(map<unsigned, unsigned>::reverse_iterator mit=mapN2iSNPMap.rbegin(); mit!=mapN2iSNPMap.rend(); ++mit)
        {
            if(mit->first==nMax)
                continue;
            ntmp=vecChips[0].sizeMapNames();
            vecChips[0].mergeEntry(vecSNPMapsX[mit->second]);
            delta+=(vecChips[0].sizeMapNames()-ntmp);
        }
        cout << "                    " << delta << " new SNPName entries added" << endl;
        
        //////// trim
        I("s", 2, "trimming maps: trimming lower density chips with HD");
        i=2;
        cout << "                    "
        << left << setfill(' ')
        << setw(6) << "chip"
        << right << setfill(' ')
        << setw(12) << "oriSize"
        << setw(12) << "newSize"
        << setw(13) << "delta"
        << endl;
        cout << "                    "
        << left << setfill(' ')
        << setw(6) << "1"
        << right << setfill(' ')
        << setw(12) << vecChips[0].sizeSetIDs()
        << setw(12) << vecChips[0].sizeSetIDs()
        << setw(13) << "0"
        << endl;
        
        for(map<unsigned, unsigned>::reverse_iterator mit=mapN2iSNPMap.rbegin(); mit!=mapN2iSNPMap.rend(); ++mit)
        {
            if(mit->first==nMax)
                continue;
            ntmp=vecSNPMapsX[mit->second].sizeSetIDs();
            vecSNPMapsX[mit->second].trimTo(vecChips[0]);
            cout << "                    "
            << left << setfill(' ')
            << setw(6) << i++
            << right << setfill(' ')
            << setw(12) << ntmp
            << setw(12) << vecSNPMapsX[mit->second].sizeSetIDs()
            << setw(13) << int(0-(ntmp-vecSNPMapsX[mit->second].sizeSetIDs()))
            << endl;
            vecChips.push_back(vecSNPMapsX[mit->second]);
        }
        //////// free memory
        vector<SNPMap>().swap(vecSNPMapsX);
        
        //////// write
        if(bOutputByChr)    // output separately by chr
        {
            I("s", 2, "writing out chromosome.data.chr separately");
            set<string>::iterator sit (vecChips[0].setSNPIDs.end());    // get the last SNPID to get the last chr
            string stmp;
            stringstream sst;
            --sit;
            stmp=(*sit).substr(0, 2);
            sst << stmp;
            sst >> maxChr;
            sst.str("");
            sst.clear();
            if(maxChr<1)
            {
                E("s", "error getting max chr number, aborting.");
                exit(-1);
            }
            char fnBuff[128];
            cout << "                   ";
            for(unsigned ichr=1; ichr<=maxChr; ++ichr)
            {
                sprintf(fnBuff, "chromosome.data.%d", ichr);
                ofstream ofs (fnBuff);
                if(!ofs)
                {
                    cout << endl;
                    E("sss", "error trying to open", fnBuff, "to write, aborting.");
                    exit(-1);
                }
                // header
                //SNPname          chrome  within  overall  location  n_chips   chip1   chip2   etc
                cout << " " << ichr << flush;
                const unsigned long n_chips(vecChips.size());
                ofs << left << setfill(' ')
                << setw(16) << "SNPname" << right
                << setw(8)  << "chrome"
                << setw(8)  << "within"
                << setw(9)  << "overall"
                << setw(10) << "location"
                << setw(9)  << "n_chips";
                for(i=0; i<n_chips; ++i)    // for each chip
                    ofs << right << setfill(' ') << setw(11) << "chip" << (i+1);
                ofs << endl;
                unsigned j, chr, pos, ilocal(1);
                vector<unsigned> vecIndexes;
                for(i=0; i<n_chips; ++i)
                    vecIndexes.push_back(1);
                i=1;

                for(sit=vecChips[0].setSNPIDs.begin(); sit!=vecChips[0].setSNPIDs.end(); ++sit) // for each SNP_ID in HD chip
                {
                    // get chr and pos
                    stmp=sit->substr(0,2);
                    sst << stmp;
                    sst >> chr;
                    sst.clear();
                    stmp=sit->substr(2);
                    sst << stmp;
                    sst >> pos;
                    sst.clear();
                    // just work with ichr==chr
                    if(ichr>chr)
                        continue;
                    else if(ichr<chr)
                        break;
                    // write the line
                    ofs << left << setfill(' ')
                    << setw(16) << (*sit) << right
                    << setw(8)  << chr
                    << setw(8)  << ilocal++
                    << setw(9)  << i++
                    << setw(10) << pos
                    << setw(9)  << n_chips;
                    for(j=0; j<n_chips; ++j)    // for each chip
                    {
                        if(vecChips[j].containsID(*sit))
                            ofs << right << setfill(' ') << setw(12) << vecIndexes[j]++;
                        else
                            ofs << right << setfill(' ') << setw(12) << 0;
                    }
                    ofs << endl;
                }
                ofs.close();
            }
            cout << endl;
        }
        else                // output all in one
        {
            I("s", 2, "writing out chromosome.data");
            ofstream ofs ("chromosome.data");
            if(!ofs)
            {
                E("s", "error trying to open 'chromosome.data' to write, aborting.");
                exit(-1);
            }
            // header
            //SNPname          chrome  within  overall  location  n_chips   chip1   chip2   etc
            const unsigned long n_chips(vecChips.size());
            ofs << left << setfill(' ')
            << setw(16) << "SNPname" << right
            << setw(8)  << "chrome"
            << setw(8)  << "within"
            << setw(9)  << "overall"
            << setw(10) << "location"
            << setw(9)  << "n_chips";
            for(i=0; i<n_chips; ++i)    // for each chip
                ofs << right << setfill(' ') << setw(11) << "chip" << (i+1);
            ofs << endl;
            set<string>::iterator sit;
            unsigned j, chr, pos, ilocal(0), lastChr(1000);
            string stmp;
            stringstream sst;
            vector<unsigned> vecIndexes;
            for(i=0; i<n_chips; ++i)
                vecIndexes.push_back(1);
            i=0.01*nMax;
            const unsigned p1p(i<2?2:i);
            i=1;
            progShow('W', 0);
            for(sit=vecChips[0].setSNPIDs.begin(); sit!=vecChips[0].setSNPIDs.end(); ++sit) // for each SNP_ID in HD chip
            {
                // get chr and pos
                stmp=sit->substr(0,2);
                sst << stmp;
                sst >> chr;
                sst.clear();
                stmp=sit->substr(2);
                sst << stmp;
                sst >> pos;
                sst.clear();
                // update ilocal
                if(lastChr!=chr)
                {
                    lastChr=chr;
                    ilocal=1;
                }
                // write the line
                ofs << left << setfill(' ')
                << setw(16) << (*sit) << right
                << setw(8)  << chr
                << setw(8)  << ilocal++
                << setw(9)  << i++
                << setw(10) << pos
                << setw(9)  << n_chips;
                for(j=0; j<n_chips; ++j)    // for each chip
                {
                    if(vecChips[j].containsID(*sit))
                        ofs << right << setfill(' ') << setw(12) << vecIndexes[j]++;
                    else
                        ofs << right << setfill(' ') << setw(12) << 0;
                }
                ofs << endl;
                if(i%p1p==0)
                    progShow('W', i/p1p);
            }
            ofs.close();
            progClear();
        }
    }
    //*********************************************************************************************
    //************************     GenConverter     ***********************************************
    //*********************************************************************************************
    GenConverter::GenConverter(vector<string> &vecFiles, PedConverter &p, MapConverter &m)
    {
        GenConverter(vecFiles, p, m, false);
    }
    
    GenConverter::GenConverter(vector<string> &vecFiles, PedConverter &p, MapConverter &m, bool bOutputByChr)
    {
        I("s", 1, "Entering GenConverter");
        unsigned i, j, k, indHighest(0), imatches;
        float mtmp, propG(0.0), propM(0.0), maxG(0.0), maxM(0.0), meanMax(0.0);
        ifstream infs;
        splitablestring aLine;
        vector<string> vecSNP_IDs;
        vector<unsigned> vecFileChipIndex;
        
        //////// matching chips
        I("s", 2, "matching chips");
        for(i=0; i<vecFiles.size(); ++i)// for each gen file, read the SNP_IDs and try matching chips
        {
            infs.open(vecFiles[i].c_str());
            if(!infs)
            {
                E("sss", "error trying to open", vecFiles[i].c_str(), "to read, aborting.");
                exit(-1);
            }
            getline(infs, aLine);// try geting lines; if it goes more than 10 lines but still not found SNP_IDs, abort
            if(aLine[0]!='#')
            {
                E("sss", "error: the given genotype file", vecFiles[i].c_str(), "is not a valid ab-genotype file: # sign not showing up in the 1st line, aborting.");
                infs.close();
                exit(-1);
            }
            vecSNP_IDs=aLine.split('#');
            if(vecSNP_IDs[0]!="SNP_IDs")
            {
                E("sss", "error: the given genotype file", vecFiles[i].c_str(), "is not a valid ab-genotype file: #SNP_IDs not showing up in the 1st line, aborting.");
                infs.close();
                exit(-1);
            }
            cout << "                    " << left << setfill(' ') << setw(36) << vecFiles[i] << right << setw(12) << vecSNP_IDs.size()-1 << " SNPs " << flush;
            meanMax=0.0;
            for(j=0; j<m.vecChips.size(); ++j)// for each chip
            {
                imatches=0;
                for(k=1; k<vecSNP_IDs.size(); ++k)  // for each SNP in the current gen file
                {
                    if(m.vecChips[j].containsName(vecSNP_IDs[k]) && m.vecChips[j].containsID(m.vecChips[j].mapName2ID[vecSNP_IDs[k]]))
                        imatches++;
                }
                //                    cout << endl << vecSNP_IDs.size()-1 << ": " << imatches << "matches / " << m.vecChips[j].sizeSetIDs() << endl;
                propG=100.0*imatches/(vecSNP_IDs.size()-1);
                propM=100.0*imatches/(m.vecChips[j].sizeSetIDs());
                mtmp=0.5*(propG+propM);
                if(mtmp>=meanMax)
                {
                    meanMax=mtmp;
                    maxM=propM;
                    maxG=propG;
                    indHighest=j;
                }
            }
            cout << "*" << right << setfill(' ') << setw(6) << setprecision(4) << maxG << "% in Chip " << (indHighest+1)
            << " [" << right << setw(9) << m.vecChips[indHighest].sizeSetIDs() << " ] * "
            << right << setfill(' ') << setw(6) << setprecision(4) << maxM << "%" << endl;
            vecFileChipIndex.push_back(indHighest);
            
            infs.close();
            vector<string>().swap(vecSNP_IDs);
        }
        
        //////// start converting
        map<string, char> mapGenCoder;
        vector<string> vecGenTmp;
        mapGenCoder["AA"]='0';
        mapGenCoder["AB"]=mapGenCoder["BA"]='1';
        mapGenCoder["BB"]='2';
        mapGenCoder["--"]=mapGenCoder["ST"]='5';
        
        set<string> setIDNotXrefed;
        map<string, unsigned> * pmapiID2nID (&p.mapiID2nID);
        ofstream ofs;
        
        if(bOutputByChr)
        {
            I("s", 2, "start converting genotype files separately by chr");
            char fnBuff[128];
            string stmp;
            stringstream sst;
            unsigned tmpChr;
            for(unsigned ichr=1; ichr<=m.maxChr; ++ichr)
            {
                sprintf(fnBuff, "genotypes.txt.UNSORTED.SORT_ME.%d", ichr);
                ofs.open(fnBuff);
                if(!ofs)
                {
                    E("sss", "error trying to open", fnBuff, "to write, aborting.");
                    exit(-1);
                }
                I("si", 3, "Chr", ichr);
                for(i=0; i<vecFiles.size(); ++i)// for each gen file, refresh all the local variables
                {
                    I("s", 4, vecFiles[i].c_str());
                    const unsigned long f1p (0.01*tellFileSize(vecFiles[i]));
                    unsigned linesRead(0);
                    infs.open(vecFiles[i].c_str());
                    SNPMap* theChip (&(m.vecChips[vecFileChipIndex[i]]));

                    map<string, unsigned> mapID2Pos;
                    unsigned idx(0);
                    set<string>::iterator sit;
                    //            I("s", 4, "filling mapID2Pos");
                    for(sit=theChip->setSNPIDs.begin(); sit!=theChip->setSNPIDs.end(); ++sit)   // fill mapID2Pos
                    {
                        stmp=(*sit).substr(0, 2);
                        sst << stmp;
                        sst >> tmpChr;
                        sst.str("");
                        sst.clear();
                        if(tmpChr==ichr)            // only output the current chr
                            mapID2Pos[*sit]=idx++;
                    }
                    map<string, unsigned> mapName2Pos;
                    map<string, string>::iterator mit;
                    //            I("s", 4, "filling mapName2Pos");
                    for(mit=theChip->mapName2ID.begin(); mit!=theChip->mapName2ID.end(); ++mit)    // fill mapName2Pos
                    {
                        if(mapID2Pos.find(mit->second)==mapID2Pos.end())    // if the name had been trimmed out
                            continue;
                        stmp=(mit->second).substr(0, 2);
                        sst << stmp;
                        sst >> tmpChr;
                        sst.str("");
                        sst.clear();
                        if(tmpChr==ichr)            // only output the current chr
                            mapName2Pos[mit->first]=mapID2Pos[mit->second];
                    }
                    
                    const unsigned nSNP (1+(unsigned)mapID2Pos.size());  // the "1" is the last char '\0'
                    char genBuffer[nSNP];

                    //            I("s", 4, "reordering SNPs and writing out");
                    getline(infs, aLine);
                    vecSNP_IDs=aLine.split('#'); // get the SNP_IDs: vecSNP_IDs=="SNP_IDs" remember!
                    const unsigned nSNPinGen(unsigned(vecSNP_IDs.size())-1);
                    vector<int> vecGenInd;
                    string gtmp;
                    int itmp;
                    //            I("s", 4, "reordering SNPs");
                    for(j=0; j<nSNPinGen+1; ++j)      // reordering SNP/genBuffer index
                    {
                        gtmp=vecSNP_IDs[j];
                        if(mapName2Pos.find(gtmp)!=mapName2Pos.end())   // this SNP should be kept - its in the map
                            vecGenInd.push_back(mapName2Pos[gtmp]);
                        else
                            vecGenInd.push_back(-1);                    // this SNP is not in the map - ignore it
                    }
                    //            I("s", 4, "for each genotype line");
                    while(getline(infs, aLine))//  for each line, do the sorting
                    {
                        if(++linesRead%150==1)
                            progShow('W', unsigned(infs.tellg()/f1p));
                        if(aLine[0]=='#')
                            continue;
                        vecGenTmp=aLine.split();
                        aLine=vecGenTmp[0];//animalID
                        if(vecGenTmp.size()!=1+nSNPinGen)
                        {
                            cout << endl;
                            E("sssisis", "line [", aLine.c_str(), "] contains", vecGenTmp.size(), "fields,", (1+nSNPinGen), "expected.");
                            continue;
                        }
                        if(pmapiID2nID->find(aLine)==pmapiID2nID->end()) //  animalID not found in the pedigree/xref
                        {
                            //                    E("sss", "animalID [", vecGenTmp[0].c_str(), "] not in the pedigree/xref, ignroing. Try puttiing it in the pedigree file next time.");
                            setIDNotXrefed.insert(aLine);
                            continue;
                        }
                        // fill genBuffer with missing '5'
                        //                I("s", 5, "filling genotype template line genBuffer");
                        for(j=0; j<nSNP-1; ++j) // the "1" is the last char '\0'
                            genBuffer[j]='5';
                        // reorder genotype and fill genBuffer
                        //                I("s", 5, "for each SNP_ID");
                        for(j=1; j<nSNPinGen+1; ++j)  // vecSNP_IDs[0] is "SNP_IDs", vecGenTmp[0] is animalID
                        {
                            if((itmp=vecGenInd[j])<0)   // ignore SNP that is not in the map
                                continue;
                            gtmp=vecGenTmp[j];
                            if(mapGenCoder.find(gtmp)!=mapGenCoder.end())
                                genBuffer[itmp]=mapGenCoder[gtmp];
                            else
                            {
                                cout << endl;
                                E("sss", "unrecognized genotype [", gtmp.c_str(), "] found, aborting.");
                                ofs.close();
                                exit(-1);
                            }
                        }
                        //                I("s", 5, "write the line");
                        genBuffer[nSNP-1]='\0';
                        // write out the current line
                        // genotypes.txt    Format: animal#   chip#     #SNPs     genotypes
                        //                          10 bytes  10 bytes  10 bytes  1 byte
                        ofs << right << setfill(' ')
                        << setw(10) << (*pmapiID2nID)[aLine]
                        << setw(10) << 1+vecFileChipIndex[i]
                        << setw(10) << nSNP-1 // no '\0' this time
                        << " " << genBuffer << endl;
                    }
                    infs.close();
                    vector<string>().swap(vecSNP_IDs);
                    progClear();
                }
                ofs.close();
            }
        }
        else
        {
            I("s", 2, "start converting genotype files");
            ofs.open("genotypes.txt.UNSORTED.SORT_ME");
            if(!ofs)
            {
                E("s", "error trying to open 'genotypes.txt.UNSORTED.SORT_ME' to write, aborting.");
                exit(-1);
            }
            for(i=0; i<vecFiles.size(); ++i)// for each gen file, refresh all the local variables
            {
                I("s", 3, vecFiles[i].c_str());
                const unsigned long f1p (0.01*tellFileSize(vecFiles[i]));
                unsigned linesRead(0);
                infs.open(vecFiles[i].c_str());
                SNPMap* theChip (&(m.vecChips[vecFileChipIndex[i]]));
                const unsigned nSNP (1+theChip->sizeSetIDs());  // the "1" is the last char '\0'
                char genBuffer[nSNP];
                map<string, unsigned> mapID2Pos;
                unsigned idx(0);
                set<string>::iterator sit;
                //            I("s", 4, "filling mapID2Pos");
                for(sit=theChip->setSNPIDs.begin(); sit!=theChip->setSNPIDs.end(); ++sit)   // fill mapID2Pos
                    mapID2Pos[*sit]=idx++;
                
                map<string, unsigned> mapName2Pos;
                map<string, string>::iterator mit;
                //            I("s", 4, "filling mapName2Pos");
                for(mit=theChip->mapName2ID.begin(); mit!=theChip->mapName2ID.end(); ++mit)    // fill mapName2Pos
                    mapName2Pos[mit->first]=mapID2Pos[mit->second];
                
                //            I("s", 4, "reordering SNPs and writing out");
                getline(infs, aLine);
                vecSNP_IDs=aLine.split('#'); // get the SNP_IDs: vecSNP_IDs=="SNP_IDs" remember!
                const unsigned nSNPinGen(unsigned(vecSNP_IDs.size())-1);
                vector<int> vecGenInd;
                string gtmp;
                int itmp;
                //            I("s", 4, "reordering SNPs");
                for(j=0; j<nSNPinGen+1; ++j)      // reordering SNP/genBuffer index
                {
                    gtmp=vecSNP_IDs[j];
                    if(mapName2Pos.find(gtmp)!=mapName2Pos.end())   // this SNP should be kept - its in the map
                        vecGenInd.push_back(mapName2Pos[gtmp]);
                    else
                        vecGenInd.push_back(-1);                    // this SNP is not in the map - ignore it
                }
                //            I("s", 4, "for each genotype line");
                while(getline(infs, aLine))//  for each line, do the sorting
                {
                    if(++linesRead%150==1)
                        progShow('W', unsigned(infs.tellg()/f1p));
                    if(aLine[0]=='#')
                        continue;
                    vecGenTmp=aLine.split();
                    aLine=vecGenTmp[0];//animalID
                    if(vecGenTmp.size()!=1+nSNPinGen)
                    {
                        cout << endl;
                        E("sssisis", "line [", aLine.c_str(), "] contains", vecGenTmp.size(), "fields,", (1+nSNPinGen), "expected.");
                        continue;
                    }
                    if(pmapiID2nID->find(aLine)==pmapiID2nID->end()) //  animalID not found in the pedigree/xref
                    {
                        //                    E("sss", "animalID [", vecGenTmp[0].c_str(), "] not in the pedigree/xref, ignroing. Try puttiing it in the pedigree file next time.");
                        setIDNotXrefed.insert(aLine);
                        continue;
                    }
                    // fill genBuffer with missing '5'
                    //                I("s", 5, "filling genotype template line genBuffer");
                    for(j=0; j<nSNP-1; ++j) // the "1" is the last char '\0'
                        genBuffer[j]='5';
                    // reorder genotype and fill genBuffer
                    //                I("s", 5, "for each SNP_ID");
                    for(j=1; j<nSNPinGen+1; ++j)  // vecSNP_IDs[0] is "SNP_IDs", vecGenTmp[0] is animalID
                    {
                        if((itmp=vecGenInd[j])<0)   // ignore SNP that is not in the map
                            continue;
                        gtmp=vecGenTmp[j];
                        if(mapGenCoder.find(gtmp)!=mapGenCoder.end())
                            genBuffer[itmp]=mapGenCoder[gtmp];
                        else
                        {
                            cout << endl;
                            E("sss", "unrecognized genotype [", gtmp.c_str(), "] found, aborting.");
                            ofs.close();
                            exit(-1);
                        }
                    }
                    //                I("s", 5, "write the line");
                    genBuffer[nSNP-1]='\0';
                    // write out the current line
                    // genotypes.txt    Format: animal#   chip#     #SNPs     genotypes
                    //                          10 bytes  10 bytes  10 bytes  1 byte
                    ofs << right << setfill(' ')
                    << setw(10) << (*pmapiID2nID)[aLine]
                    << setw(10) << 1+vecFileChipIndex[i]
                    << setw(10) << nSNP-1 // no '\0' this time
                    << " " << genBuffer << endl;
                }
                infs.close();
                vector<string>().swap(vecSNP_IDs);
                progClear();
            }
            ofs.close();
        }
        if(setIDNotXrefed.size()!=0)
        {
            E("s", "writing animalIDs [not in the pedigree file] into 'ANIMAL_IDs_NOT_IN_PEDIGREE'");
            ofs.open("ANIMAL_IDs_NOT_IN_PEDIGREE");
            for(set<string>::iterator sit=setIDNotXrefed.begin(); sit!=setIDNotXrefed.end(); ++sit)
                ofs << (*sit) << endl;
            ofs.close();
            E("s", "append animals in 'ANIMAL_IDs_NOT_IN_PEDIGREE' into your pedigree file");
        }
    }
}

int main(int argc, const char * argv[])
{
    libcbk::abg2findhap(argc, argv);
    return 0;
}
