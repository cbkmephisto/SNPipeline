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
//  VCF_Record.h
//  libHailins
//
//  Created by Hailin SU on 2/16/14.
//

#ifndef __libHailins__VCF_Record__
#define __libHailins__VCF_Record__

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include "../splitablestring.h"
#include "../appInfo.h"

using namespace std;

namespace libcbk
{
    class VCF_Record
    {
    public:
        // tab-splited
        //  0       1       2       3       4       5       6       7       8       9+
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BSHUSAM000004054458
        unsigned short  chr;
        unsigned long   pos;
        string          ID;
        vector<bool>    vec_b_GenotypeLine;

        void loadFromSS(splitablestring ss);
        void clear();
    };

    class CQueue : public vector<VCF_Record>
    {
    public:
        unsigned short curChr;
        unsigned curWindow;

        CQueue()
        {
            curChr=curWindow=0;
        }
        void clear()
        {
            // 20150603: filling the gap between two windows
            static VCF_Record prevLastRecord;
            prevLastRecord.clear();
            if(this->size())
                prevLastRecord=(*this)[this->size()-1];
            vector<VCF_Record>().swap(*this);
            curChr=curWindow=0;
            if(prevLastRecord.chr!=0)
            {
//                cout << "TEST - prevLastRecord pushed back" << endl;  // checked
                this->push_back(prevLastRecord);
            }
        }
    };
    
    class VCF_File : public appInfo
    {
    public:
        static unsigned     WINDOW_SIZE;
        const static unsigned     M1MB;
        string              fileName;
        string              breed;
        VCF_Record          lastVCFR;
        CQueue              queue;
        ifstream            infs;
        vector<string>      vecSampleIDs;
        vector<string>      sx;
        unsigned long       FILE_SIZE_1p;
        unsigned long       curTellG;

        CQueue  getNextQueue();
        CQueue  getNextQueueF(const float windowSizeMbp);
        CQueue  getNextQueueU(const unsigned nSNPperWindow);

        splitablestring     aLine;
        VCF_File(const string fn)
        {
            fileName=fn;
        }
        VCF_File(){};
        bool open();
        void close();

        void hdhdat(const string fn, const string type, const string val);
    };
}
#endif /* defined(__libHailins__VCF_Record__) */
