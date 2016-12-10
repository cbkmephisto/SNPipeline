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
//  SNPipeline_VCF2AB.h
//  SNPipeline
//
//  Created by Hailin SU on 2/26/14.
//

#ifndef SNPipeline_SNPipeline_VCF2AB_h
#define SNPipeline_SNPipeline_VCF2AB_h

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include "../splitablestring.h"
#include "../appInfo.h"

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)
#define L(...)              if(vb) lg.generalLog("I", __VA_ARGS__)

using namespace std;

namespace libcbk
{
    //****************************************************************************************** Class VCF_Record
    //****************************************************************************************** Class VCF_Record
    //****************************************************************************************** Class VCF_Record
    class VCF_Record
    {
    public:
        // tab-splited
        //  0       1       2       3       4       5       6       7       8       9+
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BSHUSAM000004054458
//        unsigned short  chr;
//        unsigned long   pos;
        string          ID;
        unsigned char   chrFlag;    // just last bit of chrString
        vector<unsigned char>   vec_c_GenotypeLine;
        vector<string>          vecStrTmp;
        unsigned length;

        VCF_Record (const unsigned l);
        void loadFromSS(splitablestring ss);
        void clear();
    };

    VCF_Record::VCF_Record(const unsigned l)
    {
        length=l;
        vec_c_GenotypeLine.resize(l);
    }

    void VCF_Record::clear()
    {
//        chr=pos=0;
//        length=0;
        ID="";
//        chrFlag=0;
        vector<unsigned char>().swap(vec_c_GenotypeLine);
        vector<string>().swap(vecStrTmp);
    }

    void VCF_Record::loadFromSS(splitablestring ss)
    {
//        vector<unsigned char>().swap(vec_c_GenotypeLine);
        static stringstream     sst;
        static unsigned         i;
        const static unsigned char uc200 (2*'0');
        static string stmp;
//        static unsigned char    ctmp;
        vecStrTmp=ss.split('\t');
        // tab-splited
        //  0       1       2       3       4       5       6       7       8       9+
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BSHUSAM000004054458
//        sst << vecStrTmp[0];
//        sst >> chr;
//        sst.clear();

//        sst << vecStrTmp[1];
//        sst >> pos;
//        sst.clear();

        chrFlag=vecStrTmp[0][vecStrTmp[0].length()-1];
        if(vecStrTmp[2][0]=='.' && vecStrTmp[2].length()==1)
        {
            ID=vecStrTmp[0];
            ID+="_";
            ID+=vecStrTmp[1];
        }
        else
            ID=vecStrTmp[2];

        for(i=9; i<vecStrTmp.size(); ++i)
        {
            // 0|1:1:0,1,0, [0]='0', [2]='1'
            //  ctmp = 0 for 0/0, 1 for 0/1 or 1/0, 2 for 1/1, 126 for ./.
/*
 vecUC2STR.push_back("./.");// 0, --
 vecUC2STR.push_back("0/0");// 1, AA
 vecUC2STR.push_back("1/1");// 2, BB
 vecUC2STR.push_back("0/1");// 3, AB
 */
//            vec_c_GenotypeLine.push_back(ctmp);
            vec_c_GenotypeLine[i-9]=vecStrTmp[i][0]+vecStrTmp[i][2]-uc200;;
        }
    }

    //****************************************************************************************** Class VCF_File
    //****************************************************************************************** Class VCF_File
    //****************************************************************************************** Class VCF_File
    class VCF_File : public appInfo
    {
    public:
        string                  fileName;
        ifstream                infs;
        vector<string>          vecSampleIDs;
        vector<string>          sx;
        vector<string>          vecSNPNames;
        vector<unsigned char>   vecGenotypes;

        unsigned long       FILE_SIZE_1p;
        bool                vb;
        bool                isValid;

        splitablestring     aLine;
        bool open(const char* fn);
        void clear();
        VCF_File(const char* fn);
    };

    bool VCF_File::open(const char* fn)
    {
        bool ret(false);
        infs.open(fn);
        if(infs)
        {
            FILE_SIZE_1p=0.01*tellFileSize(fn);
            splitablestring aLine;
            unsigned linesRead(0);      // break after reading >200 lines
            while(getline(infs, aLine)) // skip header and read sampleIDs
            {
                if(aLine.size()<2)  // make sure I can read aLine[1]
                    continue;
                else if(aLine[0]=='#' && aLine[1]=='#')  // skip meta header
                    continue;
                else if(aLine[0]=='#')  // header, read and break
                {
                    // at least 10 columns should go
                    sx=aLine.split('\t');
                    if(sx.size()<10)
                    {
                        infs.close();
                        break;              // return false
                    }
                    for(unsigned i=9; i<sx.size(); ++i)
                        vecSampleIDs.push_back(sx[i]);
                    ret=true;
                    break;
                }
                if(++linesRead>200)         // prevent large non-VCF file
                    break;
            }
        }
        fileName=fn;
        return ret;
    }

    void VCF_File::clear()
    {
        vector<string>().swap(vecSampleIDs);
        vector<string>().swap(vecSNPNames);
        vector<string>().swap(sx);
        infs.close();
    }

    VCF_File::VCF_File(const char* fn)
    {
        isValid=false;
        vb=true;
        if(!open(fn))
        {
            E("ss", "Unable to load VCF file content to proceed:", fn);
            return;
        }
        L("sss", 1, "File", fileName.c_str(), "opened successfully.");
        L("si", 2, "NumberSamples :", vecSampleIDs.size());

        VCF_Record vtmp(unsigned(vecSampleIDs.size()));

        unsigned char cf(0);

        L("s", 1, "Start loading and packing genotypes...");
        unsigned char ctmp(0), cgtmp, curCharCount(0);
        unsigned i;
        while(getline(infs, aLine)) // packing genotypes
        {
            vtmp.loadFromSS(aLine);
            vecSNPNames.push_back(vtmp.ID);


            // encoding
            // 0|1:1:0,1,0, [0]='0', [2]='1'
            //  ctmp = 0 for 0/0, 1 for 0/1 or 1/0, 2 for 1/1, 126 for ./.
            /*
             vecUC2STR.push_back("./.");// 0, --
             vecUC2STR.push_back("0/0");// 1, AA
             vecUC2STR.push_back("1/1");// 2, BB
             vecUC2STR.push_back("0/1");// 3, AB
             */
            for(i=0; i<vtmp.length; ++i)
            {
                switch(vtmp.vec_c_GenotypeLine[i])
                {
                    case 0: // 0/0, AA
                        cgtmp=1;
                        break;
                    case 1: // 0/1, AB
                        cgtmp=3;
                        break;
                    case 2: // 1/1, BB
                        cgtmp=2;
                        break;
                    default: // missing
                        cgtmp=0;
                }
                ctmp = (ctmp << 2) + cgtmp;   // left shift 2 bits and store current genotype at the lower 2 bits
                //                ctmp+=mapAB2uc[aLine];
                if(++curCharCount==4) // full, plug in and refresh
                {
                    vecGenotypes.push_back(ctmp);
                    ctmp=curCharCount=0;
                }
            }

            if(cf!=vtmp.chrFlag)
            {
                progShow('L', unsigned(infs.tellg()/FILE_SIZE_1p));
                cf=vtmp.chrFlag;
//                vtmp.clear();
//                vector<unsigned char>().swap(vtmp.vec_c_GenotypeLine);
//                vector<string>().swap(vtmp.vecStrTmp);
            }
        }

        progClear();
        if(curCharCount!=0) // somthing-left non-full
        {
            // move them to the left
            ctmp <<= 2*(4-curCharCount);
            vecGenotypes.push_back(ctmp);
        }
        L("si", 2, "NumberSNPs    :", vecSNPNames.size());

        infs.close();
        isValid=true;
    }

    //****************************************************************************************** Class SNPipeline_VCF2AB
    //****************************************************************************************** Class SNPipeline_VCF2AB
    //****************************************************************************************** Class SNPipeline_VCF2AB
    class SNPipeline_VCF2AB : public appInfo
    {
    public:
        SNPipeline_VCF2AB(const int argc, const char * argv[]);
    private:
        // methods
        void  theMain(const int argc, const char * argv[]);
        void  usage(const string& argv0);
    };
}

#endif
