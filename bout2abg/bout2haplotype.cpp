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
//  bout2haplotype.cpp
//  libHailins
//
//  Created by Hailin SU on 11/19/13.
//  Copyright (c) 2013 iastate. All rights reserved.
//
//  bout2haplotype

#include "../appInfo.h"
#include "geneticMap.h"
#include "beagleOutput.h"
#include <iomanip>

using namespace std;
namespace libcbk
{
    class bout2haplotype : public appInfo
    {
    public:
        bout2haplotype(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
        void send2write(vector<string> &vec2write, vector<string> &queue, unsigned curWin);
    };

    bout2haplotype::bout2haplotype(int argc, const char * argv[])
    {
        progName    = "bout2haplotype";
        version     = "2015-06-03 16:59 CDT";
        created     = "2013-11-19 10:59 CDT";
        updated     = "2015-06-03 16:59 CDT";
        descrip     = "read beagle output files to generate haplotype files.\n\tmatrix: 0-1 & SNP-in-col & window-grouped haplotype files";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2013-11-19 initialized
         2015-04-08 add parameter for window size
         2015-06-03 filling the gap between two adjacent windows
         */
    }

    void bout2haplotype::theMain(int argc, const char * argv[])
    {
        makeSure(argc >= 5, argv[0], "1.0 incFile mapFile beagle_(*.phased)");

        justAmap jam;
        I("s", 1, "constructing genetic map via incFile and mapFile...");
        map<short, geneticMap> gmaps = geneticMap::constructGeneticMap(argv[2], argv[3], jam);

        map<unsigned long, string>::iterator mit;
        beagleOutput bout;
        string outFileName;
        float wSize(0.0);
        
        stringstream ss;                // declaring outside the loop may cause insufficient memory
        ss << argv[1];
        ss >> wSize;
        ss.clear();
        
        unsigned long       curPOS (0);
        unsigned int        curWIN (0);
        const unsigned long WINDOW_SIZE (wSize*1000000L);  // 1M
        
        if(WINDOW_SIZE<=0)
        {
            I("sss", 1, "ERROR: window size == 1000000*[", argv[1], "]");
            return;
        }

        map<unsigned long, string> *theMap;

        for(unsigned i=4; i<argc; ++i)
        {
            I("sss", 1, "loading beagle output file : [", argv[i], "]");
            bout.clear();
            bout.load(argv[i], jam);

            ss << bout.chr;
            ss >> outFileName;
            ss.clear();
            outFileName = "bout2haplotype.chr" + outFileName;
            ofstream ofs (outFileName.c_str());
            if(!ofs)
            {
                ofs.close();
                cout << "Do not have the write access in current folder. Terminated." << endl;
                exit (-150);
            }

            vector<string> vec2write;
            vec2write.resize(2*bout.animalIDs.size());
            for(unsigned j=0; j<vec2write.size(); j++)
                vec2write[j]=bout.animalIDs[j/2];

            string              theLine, lastLine("");
            map<unsigned long, string>::iterator mit;
            vector<string> queue;
            theMap=&(gmaps[bout.chr].mapPos2Name);

            I("s", 1, "writing");
            //  go through the map (within the chr)
            curWIN=0;
            for(mit = theMap->begin(); mit != theMap->end(); ++mit)
            {
                curPOS  = mit->first;
                theLine = bout.mapMKN2Line[mit->second];
                while((1+curWIN)*WINDOW_SIZE<curPOS)
                {
                    send2write(vec2write, queue, curWIN);
                    // 20150603: fill the gap between two windows: haplotype starting from last SNP of previous window
                    if(queue.size())
                        lastLine=queue[queue.size()-1];
                    vector<string>().swap(queue);
                    curWIN++;
                    queue.push_back(lastLine);
                }
                queue.push_back(theLine);
            }
            if(queue.size()!=0)
            {
                send2write(vec2write, queue, curWIN);
                vector<string>().swap(queue);
            }
            unsigned j;
            for(j=0; j<vec2write.size(); ++j)
                ofs << vec2write[j] << endl;

            ofs.close();
            vector<string>().swap(vec2write);
        }
    }

    void bout2haplotype::send2write(vector<string> &vec2write, vector<string> &queue, unsigned int curWin)
    {
        static unsigned i, j;
        static stringstream sst;
        static string sWin, stmp;

        sst.clear();
        sst << curWin;
        sst >> sWin;
        for(j=0; j<vec2write.size(); ++j)
        {
            vec2write[j]+=" ";
            vec2write[j]+=sWin;
            vec2write[j]+="_";
        }
        // appending haplotypes from the current queue
        //nSNPs = queue.size();
        for(i=0; i<queue.size(); i++)
        {
            sst.clear();
            sst << queue[i];    // put the line of SNP haplotypes into the stream
            for(j=0; j<vec2write.size(); ++j)
            {
                sst >> stmp;
                vec2write[j] += stmp;
            }
        }
    }
}

int main(int argc, const char * argv[])
{
    libcbk::bout2haplotype(argc, argv);
    return 0;
}