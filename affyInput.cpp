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
//  affyInput.cpp
//  SNPipeline
//
//  Created by Hailin SU on 1/24/14.
//

#include "affyInput.h"
#include <cstdio>
#include <climits>
#include <iomanip>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)

using namespace std;

namespace libcbk
{
    // read basic info from affyInput file, such as
    //      numberSNPs
    //      numberSamples
    //      vecSNPIDs
    //      vecSampleIDs
    //      vecGenotypes
    // and checks for conflictions
    int affyInput::preCheck(const string fileName)
    {
        progShow('C', clock()%100);
        startTime = clock();
        //        I("s", 2, "Pre-checking...");
        infs.open(fileName.c_str());
        if(!infs)
        {
            progClear();
            E("ss", "File not found:", fileName.c_str());
            return -1;
        }

        inputFileName=fileName;
        numberSNPs = numberSamples = 0;
        string          stmp;
        splitablestring aLine;
        stringstream    sst;
//        vector<string>  vecStrx;
        unsigned        i, j;

        const unsigned long long FILE_SIZE (0.01*tellFileSize(fileName));
        // ****************** clear working space
        clear();
        vecSNPName.reserve(650000); // reserve 650k
        // ****************** read the whole file
        // 1st line: the header
        if(getline(infs, aLine))
        {
            vecStrx = aLine.split('\t');
            // vecStrx[0]==the 'probeset_id' flag
            if(vecStrx.size()/2==(1+vecStrx.size())/2)
            {
                progClear();
                E("s", "Wrong file format: column number != 1 + 2 * numberSampleIDs. Aborting.");
                infs.close();
                exit(-999);
            }

            for(i=1; i<vecStrx.size(); ++i)
            {
//                sst >> stmp;   // each sample ID appears twice
                if(vecStrx[i]!=vecStrx[i+1])
                {
                    progClear();
                    E("sss", "SampleIDs not the same for the 2 same sample columns:", vecStrx[i].c_str(), vecStrx[i+1].c_str());
                    infs.close();
                    exit(-1005);
                }
                stmp=vecStrx[++i];   // each sample ID appears twice
                if(setSampleIDs.find(stmp)==setSampleIDs.end()) // new
                {
                    vecSampleID.push_back(stmp);
                    setSampleIDs.insert(stmp);
                }
                else
                {
                    progClear();
                    E("ss", "Aborted due to duplicated Sample ID found:", stmp.c_str());
                    infs.close();
                    exit(-1006);
                }
            }

            set<string>().swap(setSampleIDs);
            numberSamples=vecSampleID.size();
            // replace whitespaces with '_'
            for(i=0; i<numberSamples; ++i)
            {
                if(vecSampleID[i].find(" ")!=string::npos)
                {
                    stmp=vecSampleID[i];
                    for(j=0; j<stmp.length(); ++j)
                        if(stmp[j]==' ') stmp[j]='_';
                    vecSampleID[i]=stmp;
                }
            }
        }
//        while (infs>>stmp)

        if(mapAB2uc.find("AB")==mapAB2uc.end()) // need init
        {
            mapAB2uc["NoCall"]=0;
            mapAB2uc["--"]=0;
            mapAB2uc["AA"]=1;
            mapAB2uc["BB"]=2;
            mapAB2uc["AB"]=3;
            mapAB2uc["BA"]=3;
        }

        unsigned curLine(1), curCharCount(0);
        unsigned char ctmp(0);
        vecGenotypes.reserve(650000*(1+numberSamples)/4);   // reserve space
        progClear();

        // rest of the file: genotype data
        while(getline(infs, aLine))
        {
            if(++curLine%6144==2) // about 1% of the 650k SNPs
                progShow('L', unsigned(infs.tellg()/FILE_SIZE));
            vecStrx = aLine.split('\t');
            if(vecStrx.size()%2==0)    // could not get the 2nd column
            {
                progClear();
                E("sssi", "Aborted due to bad data line (could not get the 2nd column for a sample) for SNP:", vecSNPName.end()->c_str(),
                  "in line", curLine);
                infs.close();
                exit(-1008);
            }
            vecSNPName.push_back(vecStrx[0]);

            for(i=1; i<vecStrx.size(); i+=2)   // vecStrx[i] is the AB genotype
            {
//                aLine=vecStrx[i];
                if(mapAB2uc.find(vecStrx[i])==mapAB2uc.end())    // not a valid genotype
                {
                    progClear();
                    E("sssi", "Aborted due to invalid SNP genotype:", vecStrx[i].c_str(), "in line", curLine);
                    infs.close();
                    exit(-1018);
                }
                /* new algorithm for storing genotype vector in the memory
                 *
                 *  AA    BB    AB/BA NoCall
                 *  1     2     3     0
                 *
                 *  each unsigned char holds up to 4 SNPs
                 *
                 *  unsigned char (5) = 0000 0101 = -- -- AA AA
                 *  unsigned char (3) = 0000 0011 = -- -- -- AB
                 *  unsigned char (?) = 1110 0100 = AB BB AA --
                 */
                ctmp = (ctmp << 2) + mapAB2uc[vecStrx[i]];   // left shift 2 bits and store current genotype at the lower 2 bits
//                ctmp+=mapAB2uc[aLine];
                if(++curCharCount==4) // full, plug in and refresh
                {
                    vecGenotypes.push_back(ctmp);
                    ctmp=curCharCount=0;
                }
            }
        }

        if(curCharCount!=0) // somthing-left non-full
        {
            // move them to the left
            //            cout << endl << hex << int(ctmp) << endl;
            ctmp <<= 2*(4-curCharCount);
            //            cout << endl << hex << int(ctmp) << endl;
            vecGenotypes.push_back(ctmp);
        }

        // check set/vec
        progClear();
        j=0.01*vecSNPName.size();
        if(j==0) j=1;
        for(i=0; i<vecSNPName.size(); ++i)
        {
            if(setSNPIDs.find(vecSNPName[i])!=setSNPIDs.end())
            {
                clear();
                E("ss", "Aborted due to duplicated SNP ID found:", vecSNPName[i].c_str());
                infs.close();
                exit(-1007);
            }
            setSNPIDs.insert(vecSNPName[i]);
            if(i%j==2) // about 1% of the SNPs
                progShow('C', unsigned(100*i/vecSNPName.size()));   // 95% max
        }
        set<string>().swap(setSNPIDs);
        numberSNPs=vecSNPName.size();
        // replace whitespaces with '_'
        for(i=0; i<numberSNPs; ++i)
        {
            if(vecSNPName[i].find(" ")!=string::npos)
            {
                stmp=vecSNPName[i];
                for(j=0; j<stmp.length(); ++j)
                    if(stmp[j]==' ') stmp[j]='_';
                vecSNPName[i]=stmp;
            }
        }
        infs.close();

        // **************** storing extra info into the object
        char buff[256];
        realpath(fileName.c_str(), buff);
        for(i=0; i<256; ++i)
        {
            if(buff[i]=='/')
                j=i;    // record the last /
            else if(buff[i]==0)
            {
                buff[j+1] = 0;  // discard fileName
                break;
            }
        }
        fullPath=buff;
        string nos;
        sst.clear();
        sst << numberSNPs;
        sst >> nos;
        sst.clear();
        ABgenotypeFileName      = "ab-genotyped-" + nos + "-Affymetrix";
        //        ABgenotypeXrefFileName  = "xref-indexed-" + content;
        ABgenotypeXrefFileName  = fullPath + "xref";
        return 1;
    }

    // basic info related to AB-genotype file
    int affyInput::postCheck()
    {
        //        I("s", 2, "Post-checking...");
        splitablestring     aLine;
        stringstream        ss;
        vector<string>      sx;
        string              stmp;
        bool                bXref(true);
        unsigned            i(0);

        progClear();
        progShow('P', 0);
        // **************** dealing with the xreFile
        infs.open(ABgenotypeXrefFileName.c_str());
        if(infs)
        {
            while(getline(infs, aLine))  // foreach line of xref
            {
                ss << aLine;            // counts the columns
                while(ss>>stmp)
                    sx.push_back(stmp);
                ss.clear();
                if(sx.size()!=2)
                {
                    progClear();
                    // wrong xref file, skipping
                    E("s", "WARNING: illegal xref file read (not a 2-column-whole-through file), skipping xref functions. #");
                    map<string, string>().swap(mapXref);
                    bXref=false;
                    break;
                }
                mapXref[sx[0]]=sx[1];
                vector<string>().swap(sx);
            }
            infs.close();
            progShow('P', 25);

            if(bXref)   // update the vecSampleID
            {
                for(i=0; i<vecSampleID.size(); ++i)
                    if(mapXref.find(vecSampleID[i])!=mapXref.end()) // got a xref
                    {
                        vecSampleID[i]=mapXref[vecSampleID[i]];
                    }
            }
        }
        progShow('P', 37);

        // **************** checks for duplication in both sample IDs and SNP IDs
        for(i=0; i<vecSampleID.size(); ++i)
            setSampleIDs.insert(vecSampleID[i]);
        progShow('P', 58);
        if(vecSampleID.size()!=setSampleIDs.size())
        {
            progClear();
            E("s", "WARNING: duplicated sampleID(s) found after xref:");
            for(i=0; i<vecSampleID.size()-1; ++i)
                for(unsigned j=i+1; j<vecSampleID.size(); ++j)
                {
                    if(vecSampleID[i]==vecSampleID[j])
                    {
                        cout << "\t\t\t" << vecSampleID[i] << endl;
                        break;
                    }
                }
            //            E("s", "Aborting... X");
            //            return -8;
        }
        set<string>().swap(setSampleIDs);
        progShow('P', 74);

        infs.open(ABgenotypeFileName.c_str());
        if(!infs)
        {
            progClear();
            I("sss", 2, "AB-genotype file [", ABgenotypeFileName.c_str(), "] not found, proceed to create A NNNNNEEEEEWWWWW ONE.");
            //            I("sis", 2, "Checking finished.", vecSampleID.size(), "new sample(s) will be added to the AB-genotype file. √");
            abFileExists = false;
            return 1;
        }
        else
        {
            abFileExists = true;
            //            I("s", 2, "AB-genotype file found, checking...");
            bool                wrongFileFormat(false);
            // **************** read the SNP_IDs: the 1st line
            if(getline(infs, aLine))
            {
                sx=aLine.split('#');
                if(sx.size()<2 || sx[0]!="SNP_IDs")
                    wrongFileFormat=true;
                // sx[1~size()-1] = SNP Name vec
            }
            else
                wrongFileFormat=true;
            if(wrongFileFormat)
            {
                progClear();
                infs.close();
                E("s", "The existing AB-genotype file has a wrong file format. Aborting. X");
                return -1;
            }
            // **************** todo: continue with the SNP_IDs
            progShow('P', 85);
            for(i=1; i<sx.size(); ++i)
            {
                if(sx[i]!=vecSNPName[i-1])
                {
                    progClear();
                    E("sssss", "ERROR: different SNP names/order found, SNP", sx[i].c_str(), "in existed ab-genotype file, while", vecSNPName[i-1].c_str(), "in the same position of affyInput file, aborting. X");
                    return -2;
                }
                //                vecABSNPName.push_back(sx[i]);
            }

            // **************** read the folder-position list in the ab-genotype file
            const unsigned long long FILE_SIZE(tellFileSize(ABgenotypeFileName.c_str())/15);
            unsigned curp(0), lastp(0);
            while(getline(infs, aLine))
            {
                if(aLine[0]!='#')
                    continue;
                vecStrx=aLine.split('#');
                if(vecStrx[0]==fullPath)
                    break;
                curp=unsigned(infs.tellg()/FILE_SIZE);
                if(curp!=lastp)
                {
                    lastp=curp;
                    progShow('P', 85+curp);
                }
            }
            if(infs.eof())   // new folder, proceed to simple appending
            {
                infs.close();
                return 1;
            }

            infs.close();
            // new sample IDs coming, doing a merging append
            return 2;
        }
    }

    // write vecGenotypes to given ofstream
    void affyInput::write(ofstream &ofs)
    {
        progClear();
        progShow('W', 0);
        static unsigned i, j, shift;
        static string stmp;
//        static unsigned char ctmp;
        if(vecUC2STR.size()==0) // need init
        {
            vecUC2STR.push_back("--");// 0
            vecUC2STR.push_back("AA");// 1
            vecUC2STR.push_back("BB");// 2
            vecUC2STR.push_back("AB");// 3
        }

        if(!abFileExists)
        {
            // write out the SNP names for the whole file
            ofs << "#SNP_IDs";
            for(i=0; i<numberSNPs; ++i)
                ofs << "#" << vecSNPName[i];
            ofs << endl;
            abFileExists=true;
        }
        // write out the full pathname of the input affyInput file
        ofs << "#" << fullPath;
        for(i=0; i<numberSamples; ++i)
            ofs << "#" << vecSampleID[i];
        ofs << endl;

        unsigned long segIndex(0);       // each char for 4 segments [ 11 01 00 10 ]
//        unsigned long charIndex(0);      // the index of vecGenotype for the current Animal & current SNP
//        unsigned char mask(0);      // the mask for current Animal for the current vecGenotypes char
        const unsigned long numx(numberSNPs*numberSamples);
        for(i=0; i<numberSamples; ++i)
        {
            if(i%7==0)  // show % every 7 samples
                progShow('W', 100*i/numberSamples);

            ofs << left << setw(32) << vecSampleID[i];
            for(j=0; j<numx; j+=numberSamples)
            {
//              segIndex =  j * numberSamples   +   i;
                //          [base]   [offset]
                segIndex =  j    +   i;
//                charIndex = segIndex / 4;           // the index of vecGenotype for the current Animal
/*
                // determine the position of genotype for the current animal
                //     2bit  3*2 0s in the left    curSeg
                shift= 2 *   (3                  - segIndex % 4);

                mask =  3 << shift;    //  0000 0011 shift left to get the mask
 
                ctmp = vecGenotypes[charIndex];
                ctmp &= mask;
                // shift ctmp right back to get the coded genotype
                ctmp >>= shift;

                // at last,
                //          mask          char = 4-seg-genotypes        shift
                //ctmp = ( (3 << shift) & vecGenotypes[charIndex] ) >> (2 * (3 - segIndex % 4));
                ofs << vecUC2STR[ctmp] << " ";

                *****************
                 *
                 *   for speed optimization
                 *   detailed algorithm see above
                 */
                ofs << vecUC2STR[( (3 << shift) & vecGenotypes[segIndex/4] ) >> (2 * (3 - segIndex % 4))] << " ";
            }
            ofs << endl;
        }
        progClear();
    }

    // affyInput file proceed to AB-genotype file with simply appending the content
    bool affyInput::AF2AB_simpleAppend()
    {
        // ****************** output
        ofstream ofs(ABgenotypeFileName.c_str(), ofstream::app | ofstream::out);  // appending mode
        if(!ofs)
        {
            progClear();
            E("s", "Could not write to file:", ABgenotypeFileName.c_str());
            return false;
        }

        write(ofs);

        ofs.close();
        I("sisssis", 2, "Finally,", numberSamples, "new sample(s) WRITEN into [", ABgenotypeFileName.c_str(), "] within",
          (clock()-startTime)/CLOCKS_PER_SEC, "seconds. √");
        return true;
    }

    // affyInput file proceed to AB-genotype file with merging appending the content
    bool affyInput::AF2AB_mergingAppend()
    {
        // creating a temporary ab-genotype file without the shared block
        //        I("s", 2, "Merging affyInput file to the existed AB-genotype file...");
        string tempFile = ABgenotypeFileName + ".TEMP.temp.tEmP.TeMp";
        ofstream opfs(tempFile.c_str());
        progClear();
        if(!opfs)
        {
            E("s", "Could not write to file:", tempFile.c_str());
            return false;
        }
        infs.open(ABgenotypeFileName.c_str());
        const unsigned long long FILE_SIZE(0.01*tellFileSize(ABgenotypeFileName));
        vector<string> sx;
        bool shouldSkip (false);
        splitablestring aLine;
        unsigned i(0), curp(0), lasp(0);
        while (getline(infs, aLine))
        {
            if(aLine[0]=='#')       // hash line, do something
            {
                sx=aLine.split('#');
                if(sx[0]==fullPath)
                    shouldSkip=true;
                else
                {
                    shouldSkip=false;
                    opfs << aLine << endl;
                }
            }
            else if(!shouldSkip)    // normal data line, just copy
                opfs << aLine << endl;
            curp=(unsigned)(infs.tellg()/FILE_SIZE);
            if(curp!=lasp)
            {
                lasp=curp;
                progShow('T', lasp);
            }
        }
        infs.close();
        progClear();

        // write the data
        write(opfs);
        opfs.close();
        
        /* rename */
        i=rename(tempFile.c_str(), ABgenotypeFileName.c_str());
        if(i==0)
        {
            I("sisssis", 2, "Finally,", numberSamples, "new sample(s) BLOCK-REPLACED into [", ABgenotypeFileName.c_str(), "] within",
                (clock()-startTime)/CLOCKS_PER_SEC, "seconds. √");
            return true;
        }
        // else
        E("ss", "ERROR: Could not rename the temp file to the AB-genotyped file:", tempFile.c_str());
        return false;
    }
    
    void affyInput::clear()
    {
        vector<string>().swap(vecSampleID);
        vector<string>().swap(vecSNPName);
        vector<string>().swap(vecStrx);
        vector<unsigned char>().swap(vecGenotypes);
        set<string>().swap(setSampleIDs);
        set<string>().swap(setSNPIDs);
        map<string, string>().swap(mapXref);
    }
}