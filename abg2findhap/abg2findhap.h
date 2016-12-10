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
//  abg2findhap.h
//  libHailins
//
//  Created by Hailin Su on 11/10/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#ifndef libHailins_abg2findhap_h
#define libHailins_abg2findhap_h

#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "../splitablestring.h"

using namespace std;

namespace libcbk
{
    //*********************************************************************************************
    //************************     abg2findhap     ************************************************
    //*********************************************************************************************
    // the main class
    class abg2findhap : public appInfo
    {
    public:
        abg2findhap(const int argc, const char * argv[]);
    private:
        bool bOutputByChr;
        void theMain(const int argc, const char * argv[]);
        void readme();
        bool parseParameters(const int argc, const char* argv[]);
        
        string pedFileName;
        vector<string> vecMapFileName;
        vector<string> vecGenFileName;
    };
    
    //*********************************************************************************************
    //************************     PedConverter     ***********************************************
    //*********************************************************************************************
    // reads pedRefiner output .csv format 3-col pedigree file, convert to findhap V3 format
    //      and make xref for IDs
    // - pedRefiner version >= 2014-11-10 so that the output is sorted
    // - date of birth will be made up
    // - unknown gender will be converted to M
    class PedConverter : public appInfo
    {
    public:
        map<string, unsigned>   mapiID2nID;
        PedConverter(const string fileName);
    };

    //*********************************************************************************************
    //************************     MapConverter     ***********************************************
    //*********************************************************************************************
    // reads all the input map files, join them, and making any xref if neccessary.
    class SNPMap
    {
    public:
        bool load(const char* fileName);
        bool containsID(const string &s);
        bool containsName(const string &s);
        unsigned sizeSetIDs();
        unsigned sizeMapNames();
        void mergeEntry(SNPMap &s); // for highest density chip, adding alias in mapName2ID
        void trimTo(SNPMap &s);     // for lower density chips, deleting not-in-common SNPs
        set<string> setSNPIDs; // output
        map<string, string> mapName2ID; // input: xref
    };
    
    class MapConverter : public appInfo
    {
    public:
        int maxChr;
        vector<SNPMap> vecChips;
        MapConverter(vector<string> &vecMapFileName, bool bOutputByChr);
        MapConverter(vector<string> &vecMapFileName);
    };
    
    //*********************************************************************************************
    //************************     GenConverter     ***********************************************
    //*********************************************************************************************
    class GenConverter : public appInfo
    {
    public:
        GenConverter(vector<string> &vecGenFileName, PedConverter &p, MapConverter &m, bool bOutputByChr);
        GenConverter(vector<string> &vecGenFileName, PedConverter &p, MapConverter &m);
    };
}
#endif
