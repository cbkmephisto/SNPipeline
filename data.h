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
//  data.h
//  SNPipeline
//
//  Created by Hailin SU on 1/29/14.
//

#ifndef __SNPipeline__data__
#define __SNPipeline__data__


#include <iostream>
#include <string>
#include "appInfo.h"
#include "splitablestring.h"
#include <iomanip>
#include <map>
#include <set>
#include <cstdlib>

using namespace std;

namespace libcbk
{
    class Data : public appInfo
    {
    public:
        /***********************************************************************
         * for both Affy and Illumina finalReport
         */
        string          strAffIll;
        string          inputFileName;
        unsigned long   numberSNPs;
        unsigned long   numberSamples;
        string          content;

        string          fullPath;
        string          ABgenotypeFileName;
        string          ABgenotypeXrefFileName;

        vector<string>          vecSNPName;
        vector<string>          vecSampleID;
        vector<unsigned char>   vecGenotypes;
        vector<string>          vecUC2STR;
        vector<string>          vecStrx;

        map<unsigned char, unsigned char>  mapAB2uc;
        map<string, string>         mapXref;

        bool                abFileExists;
        ifstream            infs;
        clock_t             startTime;

        void pipeline(const string filename);
        void checkAffIll(const string fileName);
        // checks for the AB-genotype-file
        int preCheck(const string fileName);
        int postCheck();
        // proceed to AB-genotype file
        bool DT2AB_simpleAppend();
        bool DT2AB_mergingAppend();
        void write(ofstream &ofs);

        void clear();

        /***********************************************************************
         * for Illumina finalReport
         */
        unsigned        posSNPName;
        unsigned        posSampleID;
        unsigned        posAllele1AB;
        unsigned        posAllele2AB;

        int loadI(const string fileName);
        void writeI(ofstream &ofs);

        /***********************************************************************
         * for Affymetrix input file
         */
        set<string> setSampleIDs;
        set<string> setSNPIDs;

        int loadA(const string fileName);
        void writeA(ofstream &ofs);
    };
}

#endif /* defined(__SNPipeline__data__) */
