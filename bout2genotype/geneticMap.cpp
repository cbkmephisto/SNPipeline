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
//  geneticMap.cpp
//  haplotype_diversity
//
//  Created by Hailin SU on 9/11/13.
//
#include <cstdlib>
#include "geneticMap.h"
#include "../logger.h"
#include <fstream>
using namespace std;

map<short, geneticMap> geneticMap::constructGeneticMap (string incf, string mapf, justAmap &jmp) //  chr -> geneticMap
{
    /************
                    check file access
    ************/
    cout << "            + checking file access..." << endl;
    ifstream ifinc (incf.c_str());
    if(!ifinc)
    {
        cout << "XXXXXXXX inc file not found: " << incf << endl;
        exit(-10);
    }
    ifstream ifmap (mapf.c_str());
    if(!ifmap)
    {
        cout << "XXXXXXXX map file not found: " << mapf << endl;
        ifinc.close();
        exit(-11);
    }

    map<short, geneticMap> ret;
    set<string> incSNPID;
//    unsigned long incID;
    string  incID;

    /************
                    read include file
     ************/
    cout << "            + reading include file..." << endl;
    while(ifinc >> incID)
        incSNPID.insert(incID);
    ifinc.close();
    cout << "               - " << incSNPID.size() << " reads" << endl;


    /************
                    read map file
     ************/
    cout << "            + reading map file..." << endl;
    string markerName;
//    unsigned long ISU_ID;
    unsigned short chr;
    unsigned long pos;
    geneticMap *tmap;
    stringstream ss;
    string chr_s, pos_s;
//    ifmap >> markerName >> markerName >> markerName >> markerName ; // get rid of the header line
    ifmap >> markerName >> markerName >> markerName ; // get rid of the header line
    while (ifmap)
    {
        //        ifmap>>markerName>>ISU_ID>>chr>>pos;
        ifmap>>markerName>>chr>>pos;
        // check if the current marker is in the include list
        if(incSNPID.find(markerName)==incSNPID.end()) // not found: continue
            continue;
        if(ret.find(chr)!=ret.end())    // exist chr, add to the map
        {
            tmap = &ret [chr];
            tmap->mapPos2Name.insert(pair<unsigned long, string>(pos, markerName));
        }
        else                            // new chr, create and insert
        {
            geneticMap gt;
            tmap = &gt;
            tmap->mapPos2Name.insert(pair<unsigned long, string>(pos, markerName));
            ret[chr]=*tmap;
        }
        jmp.mapMKN2chr[markerName]  = chr;
        ss << chr;
        ss >> chr_s;
        ss.clear();
        ss << pos;
        ss >> pos_s;
        ss.clear();
        jmp.mapMKN2line[markerName] = markerName + "\t" + chr_s + "\t" + pos_s;
    }
    ifmap.close();
    return ret;
}