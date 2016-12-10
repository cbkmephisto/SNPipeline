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
//  finalReport.h
//  libHailins
//
//  Created by Hailin SU on 11/12/13.
//
//  process finalReport files to write/append to the appropriate files.
//  output files were named after the Content line of the finalReport file.

#ifndef __libHailins__finalReport__
#define __libHailins__finalReport__

#include <iostream>
#include <string>
#include "appInfo.h"
#include "splitablestring.h"
#include <iomanip>
#include <map>
#include <set>
#include <stdlib.h>

using namespace std;

namespace libcbk
{
    class finalReport : public appInfo
    {
    public:
        // basic info from finalReport file filled in preCheck
        string          inputFileName;
        unsigned        numberSNPs;
        unsigned        numberSamples;
        string          content;

        unsigned        posSNPName;
        unsigned        posSampleID;
        unsigned        posAllele1AB;
        unsigned        posAllele2AB;

        string          finalReportFullPath;
        string          ABgenotypeFileName;
        string          ABgenotypeXrefFileName;

        vector<string>  vecSNPName;
        vector<string>  vecSampleID;
        vector<string>  vecStrx;

        // basic info from related AB-genotype file filled in postCheck
        bool                abFileExists;
//        vector<string>      vecABSNPName;
//        vector<string>      vecSampleID2bAppended;
//        set<string>         setABSampleID;
        map<string, string> mapXref;

        clock_t startTime;
        map<string, unsigned char> mapAB2uc;
        vector<string> vecUC2STR;

        // read basic info from finalReport file, such as
        //      numberSNPs
        //      numberSamples
        //      content
        //      vecSNPIDs
        //      vecSampleIDs
        //      posSNPName
        //      posSampleID
        //      posAllele1AB
        //      posAllele2AB
        int preCheck(const string fileName);
        // and checks for the AB-genotype-file
        int postCheck();
        // finalReport file proceed to AB-genotype file
        bool FR2AB_simpleAppend();
        bool FR2AB_mergingAppend();

//        bool proc(const string fileName);

//        void print();
//        void print(ofstream ofs);
        void clear();
    };
}

#endif /* defined(__libHailins__finalReport__) */
