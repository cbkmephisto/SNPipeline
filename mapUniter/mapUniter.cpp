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
//  mapUniter.cpp
//  libHailins
//
//  Created by Hailin SU on 3/24/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#include "mapUniter.h"

//  mapUniter

using namespace std;
namespace libcbk
{
    class mapUniter : public appInfo
    {
    public:
        mapUniter(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
    };

    mapUniter::mapUniter(int argc, const char * argv[])
    {
        progName    = "mapUniter";
        version     = "2014-03-24 17:32 CDT";
        created     = "2014-03-24";
        updated     = "N/A";
        descrip     = "Unite map files to get a snp_info file for FImpute.";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2013-03-24 init
         */
    }

    void mapUniter::theMain(int argc, const char * argv[])
    {
        makeSure(argc >= 3, argv[0], "map1 map2 (*)");

        map<unsigned long, vector<SNP> > mapOfMaps;  // map < nSNPs - vector >
        unsigned i, j, p1;

        for(unsigned argind=1; argind<argc; ++argind)
        {
            I("sss", 1, "{", argv[argind], "}");
            ifstream infs(argv[argind]);
            if(!infs)
            {
                E("ss", "Unable to open file", argv[argind]);
                continue;
            }
            vector<SNP> vecSNP;

            // **************** reading the map File
            splitablestring aLine;
            string stmp;
            stringstream ss;
            unsigned long ltmp;
            unsigned utmp;
            vector<string> sx;
            unsigned long FILE_SIZE(0.01*tellFileSize(argv[argind]));
            if(FILE_SIZE<100)
                FILE_SIZE=100;
            unsigned linesRead(0);
            set<string> setSNPName;
            while(getline(infs, aLine))  // foreach line of mapFile
            {
                sx=aLine.split();
                if(sx.size()!=3)
                {
                    // wrong map file
                    cout << endl;
                    E("s", "ERROR: illegal map file read (not a 3-column-whole-through space/tab-delimited file). Skipping. X");
                    cerr << "              line content: " << aLine << endl;
                    infs.close();
                    break;
                }
                // skip # line, for header
                if(sx[0][0]=='#')
                    continue;

                ss.str("");
                ss << sx[1];
                ss >> utmp;     // chr
                ss.clear();
                ss << sx[2];
                ss >> ltmp;     // pos
                ss.clear();

                // skip chr0, pos0 pool SNP
                if(utmp==0 || 0==ltmp)
                    continue;
                
                // duplicated SNP_ID: only allow 1st records
                if(setSNPName.find(sx[0])==setSNPName.end())
                {
                    vecSNP.push_back(SNP(sx[0], utmp, ltmp));
                    setSNPName.insert(sx[0]);
                }
                if(++linesRead%1000==1)
                    progShow('L', unsigned(infs.tellg()/FILE_SIZE));
            }
            set<string>().swap(setSNPName);
            infs.close();   // finish reading the mapFile
            sort(vecSNP.begin(), vecSNP.end(), SNP::SNPcomp);
            for(i=0; i<vecSNP.size()-1; ++i)
                if(vecSNP[i].equalsbyPos(vecSNP[i+1]))
                {
                    ++vecSNP[i+1].pos;
                    if(i<vecSNP.size()-2 && vecSNP[i+1].equalsbyPos(vecSNP[i+2]))
                        ++vecSNP[i+2].pos;
                }
            mapOfMaps[vecSNP.size()]=vecSNP;
            progClear();
            I("is", 2, vecSNP.size(), "SNPs read from the map.");
        }

        // **************** write all the maps out
        // SNP_ID           Chr      Pos   Chip1   Chip2
        ofstream ofs ("snp_info.txt");
        if(!ofs)
        {
            E("s", "Initializing output file { snp_info.txt } failed. Aborting. X");
            return;
        }

        // header
        ofs << setfill(' ') << left
        << setw(55) << "SNP_ID" << right
        << setw(5)  << "Chr"
        << setw(12) << "Pos";
        for(i=0; i<mapOfMaps.size(); ++i)
            ofs << setfill(' ') << right << setw(13) << "Chip" << i+1;
        ofs << endl;

        vector < map<unsigned long, vector<SNP> >::reverse_iterator > vecPtSNP; // vecPtSNP[i]->second is the vector of SNP
        vector < unsigned > vecInd, vecLocalInd;
        map<unsigned long, vector<SNP> >::reverse_iterator xit;

        for(xit=mapOfMaps.rbegin(); xit!=mapOfMaps.rend(); ++xit)
        {
            vecPtSNP.push_back( xit );
            vecInd.push_back(0);            // for output
            vecLocalInd.push_back(0);       // for input
        }

        SNP *si, *sj;
        p1=unsigned(vecPtSNP[0]->second.size()*0.01);
        if(p1<1) p1=10;
        for(i=0; i<vecPtSNP[0]->second.size(); ++i) // for each SNP in the HD map
        {
            si = &(vecPtSNP[0]->second[i]);
            if(si->chr==31)     // skip chr31
                continue;
            ofs << setfill(' ') << left
            << setw(55) << si->ID << right
            << setw(5)  << si->chr
            << setw(12) << si->pos
            << setw(14) << ++vecInd[0];
            for(j=1; j<vecInd.size(); ++j)  // for each lower density chip
            {
                sj=&(vecPtSNP[j]->second[vecLocalInd[j]]);   // sj is the pointer pointing to the next SNP of current chip
                
                // locate sj
                // keep moving if
                //  - sj is before si &&
                //  - vecInd[j] in safe range &&
                //  - sj is not si

                while(si->gt(*sj) && vecInd[j]<=vecPtSNP[j]->second.size()-1 && !si->equalsbyID(*sj))
                    sj=&(vecPtSNP[j]->second[++vecLocalInd[j]]);

                if(si->equalsbyID(*sj))    // matches
                    ofs << setfill(' ') << right << setw(14) << ++vecInd[j];
//                else if(si->equalsbyPos(*sj))
//                    ofs << setfill(' ') << right << setw(14) << ++vecInd[j];
                else
                    ofs << setfill(' ') << right << setw(14) << 0;
            }
            ofs << endl;
            if(i%p1==0)
                progShow('W', i/p1);
        }
        progClear();
        ofs.close();

    }
}

int main(int argc, const char * argv[])
{
    libcbk::mapUniter(argc, argv);
    return 0;
}