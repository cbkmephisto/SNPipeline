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
//  TrimDown2Map.cpp
//  libHailins
//
//  Created by Hailin SU on 10/3/13.
//

//  TrimDown2Map
#include "../appInfo.h"
#include "../splitablestring.h"
#include <set>
#include <vector>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)

using namespace std;
namespace libcbk
{
    class TrimDown2Map : public appInfo
    {
    public:
        TrimDown2Map(const int argc, const char * argv[]);
    private:
        void theMain(const int argc, const char * argv[]);
    };

    TrimDown2Map::TrimDown2Map(const int argc, const char * argv[])
    {
        progName    = "TrimDown2Map";
        version     = "2014-01-07 13:16 CDT";
        created     = "2014-01-07";
        updated     = "N/A";
        descrip     = "TrimDown2Map - trim the given genotype file down to a mapped one.";
        //        displayInfo();
        //        lg.initLogger();
        theMain(argc, argv);
        //        lg.finalizeLogger();
        /*
         updates
         2014-01-07 init
         */
    }

    void TrimDown2Map::theMain(const int argc, const char * argv[])
    {
        makeSure(argc == 3, argv[0], "mapFile genotypeFile");
        ifstream infs (argv[2]);
        if(!infs)
        {
            E("s", "wrong genotype file");
            exit(0);
        }
        infs.close();
        ofstream ofs((string("tdmout.")+argv[1]+"."+argv[2]).c_str());
        if(!ofs)
        {
            E("s", "no write access");
            exit(0);
        }
        infs.open(argv[1]);
        if(!infs)
        {
            ofs.close();
            E("s", "wrong map file");
            exit(0);
        }

//        I("s", 1, "reading the map");
        // read the map
        splitablestring aLine;
        set<string> setValideSNPNames;

        unsigned short chr;
        unsigned long  pos;
        stringstream sst;
        while(getline(infs, aLine))// && aLine[0]=='#')
        {
            sst << aLine;
            sst >> aLine >> chr >> pos;   // get the 1st col
            sst.str("");
            sst.clear();
            if(chr!=0 && pos!=0)
                setValideSNPNames.insert(aLine);
        }
//        else
/*        {
            infs.close();
            ofs.close();
            E("s", "Wrong map file content. The 1st line should be the header line and starts with '#'.");
            exit(0);
        }
*/
        infs.close();
        // setValidateSNPNames ready.

        infs.open(argv[2]);// open genotype file
        // init the indexed result
        vector<unsigned long> vecIndexShouldBeKept;
        vector<string> vecLine;
        getline(infs, aLine);// get the header line
        vecLine=aLine.split(" #");
//        vecShouldBeKept.resize(vecLine.size());
        unsigned i;
//        vecShouldBeKept[0]=true;    // 1st column: keep always
        vecIndexShouldBeKept.push_back(0);
        // fills the vecShouldBeKept
        for(i=1; i<vecLine.size(); ++i)
            if(setValideSNPNames.find(vecLine[i])!=setValideSNPNames.end()) // current column SNP NAME in map, keep it
                vecIndexShouldBeKept.push_back(i);

        // write out the header line
        for(i=0; i<vecIndexShouldBeKept.size(); ++i)
            ofs << "#" << vecLine[vecIndexShouldBeKept[i]];
        ofs << endl;

        // write out the whole file
//        I("sss", 1, "{", argv[2], "}");
        unsigned long FILE_SIZE (0.01*tellFileSize(argv[2]));
        if(FILE_SIZE<100)
            FILE_SIZE=100;
        unsigned linesRead(0);
        string oLine;
        while(getline(infs, aLine)) // for each line
        {
            if(aLine[0]=='#')
                continue;           // 20141113: skip # lines
            vecLine=aLine.split(' ');   // split
            oLine=vecLine[vecIndexShouldBeKept[0]];
            for(i=1; i<vecIndexShouldBeKept.size(); ++i)
            {
                oLine+=" ";
                oLine+=vecLine[vecIndexShouldBeKept[i]];
            }
            ofs << oLine << endl;
            if(++linesRead%100==1)
            {
                vector<string>().swap(vecLine);
                progShow('T', unsigned(infs.tellg()/FILE_SIZE));
            }
        }
        progClear();
        infs.close();
        ofs.close();
//        I("s",1,"done");
    }
}

int main(int argc, const char * argv[])
{
    libcbk::TrimDown2Map(argc, argv);
    return 0;
}
