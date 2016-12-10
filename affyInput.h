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
//  affyInput.h
//  SNPipeline
//
//  Created by Hailin SU on 1/24/14.
//

#ifndef __SNPipeline__affyInput__
#define __SNPipeline__affyInput__

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
    class affyInput : public appInfo
    {
    public:
        // basic info from affyInput file filled in preCheck
        string          inputFileName;
        unsigned long   numberSNPs;
        unsigned long   numberSamples;

        string          fullPath;
        string          ABgenotypeFileName;
        string          ABgenotypeXrefFileName;

        vector<string>          vecSNPName;
        vector<string>          vecSampleID;
//        vector<string>  vecGenotypes;
        vector<unsigned char>   vecGenotypes;

        set<string> setSampleIDs;
        set<string> setSNPIDs;

        map<string, unsigned char> mapAB2uc;
        vector<string> vecUC2STR;
        vector<string>  vecStrx;
        // basic info from related AB-genotype file filled in postCheck
        bool                abFileExists;
        map<string, string> mapXref;

        clock_t startTime;

        ifstream        infs;
        // read basic info from finalReport file, such as
        //      numberSNPs
        //      numberSamples
        //      vecSNPIDs
        //      vecSampleIDs
        int preCheck(const string fileName);
        // and checks for the AB-genotype-file
        int postCheck();
        // finalReport file proceed to AB-genotype file
        bool AF2AB_simpleAppend();
        bool AF2AB_mergingAppend();
        void load();
        void write(ofstream &ofs);
        void clear();
    };
}

#endif /* defined(__SNPipeline__affyInput__) */
