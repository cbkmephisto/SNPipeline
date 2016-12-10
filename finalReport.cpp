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
//  finalReport.cpp
//  libHailins
//
//  Created by Hailin SU on 11/12/13.
//

#include "finalReport.h"
#include <cstdio>
#include <climits>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)

using namespace std;

namespace libcbk
{
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
    // and checks for conflictions
    int finalReport::preCheck(const string fileName)
    {
        startTime = clock();
//        I("s", 2, "Pre-checking...");
        ifstream infs (fileName.c_str());
        if(!infs)
        {
            E("ss", "File not found:", fileName.c_str());
            return -1;
        }

        inputFileName=fileName;
        numberSNPs = numberSamples = UINT_MAX;
        content = "___";

        splitablestring dataLine;
        unsigned i (0), j;
        string aLine;

        // **************** [Header] part
        while(infs >> aLine && aLine!="[Data]")
        {
//            if(aLine=="[Header]")
//                I("s", 2, aLine.c_str());
            if(++i>200) break;       // avoid large file
            if(aLine=="Content")
            {
                infs >> content;
//                I("ss", 3, "Content       = ", content.c_str());
            }
            else if(aLine == "Num")
            {
                infs >> aLine;
                if(aLine == "SNPs")
                {
                    infs >> numberSNPs;
//                    I("si", 3, "numberSNPs    = ", numberSNPs);
                }
                else if(aLine == "Samples")
                {
                    infs >> numberSamples;
//                    I("si", 3, "numberSamples = ", numberSamples);
                }
            }
        }
        if(numberSNPs==UINT_MAX || numberSamples==UINT_MAX || content=="___")
        {
            infs.close();
            E("s", "ERROR: wrong header part of the header section, maybe the specified file is NOT a finalReport file. Aborting. X");
            return -2;
        }
        if(aLine!="[Data]")
        {
            infs.close();
            E("s", "ERROR: no [Data] section found, maybe the specified file is NOT a finalReport file. Aborting. X");
            return -3;
        }
//        I("s", 2, aLine.c_str());
        getline(infs, aLine);   // the \n after [Data] line

        // **************** header of [Data] part
        getline(infs, dataLine);
        vecStrx = dataLine.split('\t');

        posSNPName = posSampleID = posAllele1AB = posAllele2AB = UINT_MAX;

        for(i=0; i<vecStrx.size(); ++i)
        {
            if(vecStrx[i]=="SNP Name")
            {
                posSNPName   = i;
//                I("si", 3, "posSNPName    = ", posSNPName);
            }
            else if(vecStrx[i]=="Sample ID")
            {
                posSampleID  = i;
//                I("si", 3, "posSampleID   = ", posSampleID);
            }
            else if(vecStrx[i]=="Allele1 - AB")
            {
                posAllele1AB = i;
//                I("si", 3, "posAllele1AB  = ", posAllele1AB);
            }
            else if(vecStrx[i]=="Allele2 - AB")
            {
                posAllele2AB = i;
//                I("si", 3, "posAllele2AB  = ", posAllele2AB);
            }
        }
        if(posSNPName==UINT_MAX || posSampleID==UINT_MAX || posAllele1AB==UINT_MAX || posAllele2AB==UINT_MAX)
        {
            infs.close();
            E("s", "ERROR: wrong header structure of the [Data] section. Aborting... X");
            return -4;
        }

        // **************** going through the data section
        if(mapAB2uc.find("AB")==mapAB2uc.end()) // need init
        {
            mapAB2uc["NoCall"]=0;
            mapAB2uc["--"]=0;
            mapAB2uc["AA"]=1;
            mapAB2uc["BB"]=2;
            mapAB2uc["AB"]=3;
            mapAB2uc["BA"]=3;
        }
        // reads all the raw-Sample/SNP IDs
        vector<string> vecRawSampleIDs;
        vector<string> vecRawSNPIDs;
        const unsigned long long FILE_SIZE (tellFileSize(fileName));
        const unsigned vecSize (numberSamples*numberSNPs);
        const unsigned p1p     (0.015*vecSize);
        vecRawSampleIDs.reserve(vecSize);
        vecRawSNPIDs.reserve(vecSize);
        unsigned lastp(0);
        unsigned char ctmp(0);

        while(getline(infs, dataLine))
        {
            if(++lastp%p1p==0)
                progShow('l', (unsigned)(100*infs.tellg()/FILE_SIZE));
            vecStrx=dataLine.split('\t');
            vecRawSNPIDs.push_back(vecStrx[0]);
            vecRawSampleIDs.push_back(vecStrx[1]);
            // encoding
        }
        progClear();
        infs.close();

        // check raw-ids
        string prevID("--"), sidtmp;
        unsigned curIDindex(0), SNPIDx(0);
        for(i=0; i<vecSize; ++i)
        {
            if(prevID!=vecRawSampleIDs[i]) // new Animal
            {
                prevID=vecRawSampleIDs[i];
                progShow('C', 100*curIDindex/numberSamples);
                ++curIDindex;

                if(curIDindex>1)// when SNP name reading has been done,
                {
                    SNPIDx=0;
                    if(curIDindex==2)   // collect the info
                    {
                        if(vecSNPName.size()!=numberSNPs)
                        {
                            progClear();
                            I("sis", 3, "SNP IDs read  = ", vecSNPName.size(), "X");
                            E("s", "ERROR: SNP_IDs_read != numberSNPs from [Header] section, bad finalReport file, aborting. X");
                            return -5;
                        }
                    }
                }
                vecSampleID.push_back(prevID);
            }
            if(curIDindex==1)   // data of 1st sample: fill SNP Names
            {
                vecSNPName.push_back(vecRawSNPIDs[i]);
            }
            else                // data of others: check the SNP IDs
            {
                if(vecSNPName[SNPIDx++]==vecRawSNPIDs[i])
                    continue;
                progClear();
                E("sssssss", "ERROR: different SNP names/order found, Sample", prevID.c_str(),
                  "- SNP", vecRawSNPIDs[i].c_str(), "should be", vecSNPName[SNPIDx-1].c_str(), ", bad finalReport file, aborting. X");
                return -6;
            }
        }

        vector<string>().swap(vecRawSNPIDs);
        vector<string>().swap(vecRawSampleIDs);

        if(vecSampleID.size()!=numberSamples)
        {
            progClear();
            I("sis", 3, "number of SampleIDs read = ", vecSampleID.size(), "X");
            E("s", "ERROR: Sample_IDs_read != numberSamples from [Header] section, bad finalReport file, aborting. X");
            return -7;
        }
        if(vecSNPName.size()!=numberSNPs)
        {
            progClear();
            I("sis", 3, "number of SNPNames read = ", vecSNPName.size(), "X");
            E("s", "ERROR: SNP_Names_read != numberSNPs from [Header] section, bad finalReport file, aborting. X");
            return -7;
        }
        // replace whitespaces with '_'
        string stmp;
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

        // **************** checks for duplication in both sample IDs and SNP IDs
        set<string> setSampleIDs, setSNPNames;
        for(i=0; i<vecSampleID.size(); ++i)
            setSampleIDs.insert(vecSampleID[i]);
        if(vecSampleID.size()!=setSampleIDs.size())
        {
            progClear();
            E("s", "Duplicated sampleID(s) found:");
            for(i=0; i<vecSampleID.size()-1; ++i)
                for(unsigned j=i+1; j<vecSampleID.size(); ++j)
                {
                    if(vecSampleID[i]==vecSampleID[j])
                    {
                        cout << "\t\t\t" << vecSampleID[i] << endl;
                        break;
                    }
                }
            E("s", "Aborting... X");
            return -8;
        }
        set<string>().swap(setSampleIDs);

        for(i=0; i<vecSNPName.size(); ++i)
            setSNPNames.insert(vecSNPName[i]);
        if(vecSNPName.size()!=setSNPNames.size())
        {
            progClear();
            E("s", 3, "Duplicated SNP Name(s) found in the finalReport file.");
            for(i=0; i<vecSNPName.size()-1; ++i)
                for(unsigned j=i+1; j<vecSNPName.size(); ++j)
                {
                    if(vecSNPName[i]==vecSNPName[j])
                    {
                        cout << "\t\t\t" << vecSNPName[i] << endl;
                        break;
                    }
                }
            E("s", 3, "Aborting... X");
            return -9;
        }
        set<string>().swap(setSNPNames);
//        I("s", 3, "This finalReport file has been validated to be legal. √");

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
        finalReportFullPath=buff;
/*        if(finalReportFullPathName.find(" ")!=string::npos || finalReportFullPathName.find("\t")!=string::npos)
            for(i=0; i<finalReportFullPathName.length(); ++i)
                if(finalReportFullPathName[i]==' ' || finalReportFullPathName[i]=='\t')
                    finalReportFullPathName[i]='_';
 */
        stringstream ss;
        string nos;
        ss << numberSNPs;
        ss >> nos;
        ss.clear();
        ABgenotypeFileName      = "ab-genotyped-" + nos + "-" + content;
//        ABgenotypeXrefFileName  = "xref-indexed-" + content;
        ABgenotypeXrefFileName  = finalReportFullPath + "xref";
        progClear();
        return 1;
    }

    // basic info related to AB-genotype file
    int finalReport::postCheck()
    {
//        I("s", 2, "Post-checking...");
        splitablestring     aLine;
        stringstream        ss;
        vector<string>      sx;
        string              stmp;
        bool                bXref(true);
        unsigned            i(0);

        // **************** dealing with the xreFile
        ifstream ifs(ABgenotypeXrefFileName.c_str());
        const unsigned long long FILE_SIZE(0.01*tellFileSize(ABgenotypeXrefFileName));
        if(ifs)
        {
            while(getline(ifs, aLine))  // foreach line of xref
            {
                progShow('P', (unsigned)(ifs.tellg()/FILE_SIZE));
                ss << aLine;            // counts the columns
                while(ss>>stmp)
                    sx.push_back(stmp);
                ss.clear();
                if(sx.size()!=2)
                {
                    // wrong xref file, skipping
                    E("s", "WARNING: illegal xref file read (not a 2-column-whole-through file), skipping xref functions. #");
                    map<string, string>().swap(mapXref);
                    bXref=false;
                    break;
                }
                mapXref[sx[0]]=sx[1];
                vector<string>().swap(sx);
            }
            ifs.close();
            progClear();
            if(bXref)   // update the vecSampleID
            {
                for(; i<vecSampleID.size(); ++i)
                    if(mapXref.find(vecSampleID[i])!=mapXref.end()) // got a xref
                    {
                        vecSampleID[i]=mapXref[vecSampleID[i]];
                    }
            }
        }

        // **************** checks for duplication in both sample IDs and SNP IDs
        set<string> setSampleIDs;
        for(i=0; i<vecSampleID.size(); ++i)
            setSampleIDs.insert(vecSampleID[i]);
        if(vecSampleID.size()!=setSampleIDs.size())
        {
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

        ifs.open(ABgenotypeFileName.c_str());
        if(!ifs)
        {
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
            if(getline(ifs, aLine))
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
                ifs.close();
                E("s", "The existing AB-genotype file has a wrong file format. Aborting. X");
                return -1;
            }
            // **************** todo: continue with the SNP_IDs
            for(i=1; i<sx.size(); ++i)
            {
                if(sx[i]!=vecSNPName[i-1])
                {
                    E("sssss", "ERROR: different SNP names/order found, SNP", sx[i].c_str(), "should be", vecSNPName[i-1].c_str(), ", bad finalReport file, aborting. X");
                    return -2;
                }
//                vecABSNPName.push_back(sx[i]);
            }

            // **************** read the folder-position list in the ab-genotype file
            while(getline(ifs, aLine, '#') && aLine!=finalReportFullPath);
            if(ifs.eof())   // new folder, proceed to simple appending
            {
                ifs.close();
                return 1;
            }

            // **************** read the Animal IDs
//            getline(ifs, aLine);    // the rest of the line contains all the animal IDs
            //vector<string> vecABSampleID = aLine.split('#');
//            sx=aLine.split('#');
//            for(i=0; i<sx.size(); ++i)
//                setABSampleID.insert(sx[i]);    // setABSampleID just includes the animal IDs in the hash line with the same folder of the current processing finalReport file
/*          
 *
 *          NO SKIPPING... EVERYTHING WILL BE RE-PROCESSED (RETURN 2) AFTER
 *                  DEC 17, 2013
 *
 *
            // **************** checking sample IDs in finalReport with AB-genotype file
            for(i=0; i<vecSampleID.size(); ++i)
            {
                // store the sampleIDs that only appeared in the finalReport file
                // shared sampleIDs from the new finalReport file will be discard !!!!!!!
                if(setABSampleID.find(vecSampleID[i])==setABSampleID.end())
                    vecSampleID2bAppended.push_back(vecSampleID[i]);
            }
            // if all the sample IDs in the finalReport file are in the AB-genotype file, abort
            if(vecSampleID2bAppended.size()==0)
            {
                ifs.close();
                I("s", 2, "AB-genotype file contains all the sample IDs in the finalReport file. Skipping. *");
                return 0;
            }
*/
//          I("sisis", 2, "Checking finished.", j, "/", vecSampleID.size(), "new sample(s) will be added to the AB-genotype file. √");
            ifs.close();
            // new sample IDs coming, doing a merging append
            return 2;
        }
    }

    // finalReport file proceed to AB-genotype file with simply appending the content
    // write out the AB genotype file while reading finalReport file
    bool finalReport::FR2AB_simpleAppend()
    {
        ifstream infs (inputFileName.c_str());
        splitablestring dataLine;
        vector<string> sx;
        unsigned i (0);

        string aLine, Content ("___");

        while(infs >> aLine && aLine!="[Data]");
        getline(infs, aLine);   // the \n after [Data] line
        getline(infs, aLine);   // the header line of [Data] section

//        I("s", 2, "Writing to the output AB-genotype file...");
        // going through the data section
        string prevID ("--"), curAnimID;
        unsigned curIDindex(0);
        string lineGenotype;
        lineGenotype.reserve(numberSNPs*3);     // reserve space for efficiency

        ofstream opfs(ABgenotypeFileName.c_str(), ofstream::app | ofstream::out);  // appending mode
        if(!opfs)
        {
            infs.close();
            E("s", "Could not write to file:", ABgenotypeFileName.c_str());
            return false;
        }

        // write out the SNP names for the whole file
        if(!abFileExists)
        {
            opfs << "#SNP_IDs";
            for(i=0; i<vecSNPName.size(); ++i)
                opfs << "#" << vecSNPName[i];
            opfs << endl;
            abFileExists=true;
        }

        // write out the full pathname of the input finalReport file
        opfs << "#" << finalReportFullPath;
        for(i=0; i<vecSampleID.size(); ++i)
            opfs << "#" << vecSampleID[i];
        opfs << endl;

        while(getline(infs, dataLine))
        {
            sx=dataLine.split('\t');
            if(prevID!=sx[1]) // new Animal
            {
                prevID=sx[1];
                progShow('W', 100*++curIDindex/numberSamples);

                if(curIDindex>1)// write
                    opfs
                    << left << setw(32) << curAnimID
                    << lineGenotype << endl;
                // renew vars
                lineGenotype="";
                curAnimID=sx[1];
                // replace whitespaces with '_'
                if(curAnimID.find(" ")!=string::npos)
                    for(i=0; i<curAnimID.length(); ++i)
                        if(curAnimID[i]==' ') curAnimID[i]='_';
                if(mapXref.find(curAnimID)!=mapXref.end())
                    curAnimID=mapXref[curAnimID];
            }
            lineGenotype+=sx[posAllele1AB];
            lineGenotype+=sx[posAllele2AB];
            lineGenotype+=" ";
        }
        // write out the data of the last animal
        opfs
        << left << setw(32) << curAnimID
        << lineGenotype << endl;

        infs.close();
        opfs.close();
        progClear();
        I("sisssis", 2, "Finally,", vecSampleID.size(), "new sample(s) WRITEN into [", ABgenotypeFileName.c_str(), "] within",
          (clock()-startTime)/CLOCKS_PER_SEC, "seconds. √");
        return true;
    }

    // finalReport file proceed to AB-genotype file with merging appending the content
    bool finalReport::FR2AB_mergingAppend()
    {
        // creating a temporary ab-genotype file without the shared block
//        I("s", 2, "Merging finalReport file to the existed AB-genotype file...");
        string tempFile = ABgenotypeFileName + ".TEMP.temp.tEmP.TeMp";
        ofstream opfs(tempFile.c_str());
        if(!opfs)
        {
            E("s", "Could not write to file:", tempFile.c_str());
            return false;
        }
        ifstream infs (ABgenotypeFileName.c_str());
//        vector<string> vecReservedBlock;
        vector<string> sx;
        bool shouldSkip (false);
        splitablestring aLine;
        unsigned i(0);
        const unsigned long long FILE_SIZE(0.01*tellFileSize(ABgenotypeFileName));
        while (getline(infs, aLine))
        {
            progShow('T', (unsigned)(infs.tellg()/FILE_SIZE));
            if(aLine[0]=='#')       // hash line, do something
            {
                sx=aLine.split('#');
                if(sx[0]==finalReportFullPath)
                {
                    shouldSkip=true;
                    // reconstruct the hash line later, do not reserve original hash line.
//                    for(i=0; i<vecSampleID2bAppended.size(); ++i)
//                    {
//                        aLine+="#";
//                        aLine+=vecSampleID2bAppended[i];
//                    }
//                    vecReservedBlock.push_back(aLine);
                }
                else
                {
                    shouldSkip=false;
                    opfs << aLine << endl;
                }
            }
            else if(!shouldSkip)    // normal data line, just copy
                opfs << aLine << endl;
/*          JUST SKIP!!!! AFTER DEC 17, 2013

            else                    // push me to the reserved block only if I'm valid (appears in the 2nd column of xref file)
            {
                sx=aLine.split(' ');
                if(setValidOutputABSampleID.find(sx[0])!=setValidOutputABSampleID.end())
                    vecReservedBlock.push_back(aLine);
            }
*/
        }

        infs.close();
        infs.open(inputFileName.c_str());
        splitablestring dataLine;

        while(infs >> aLine && aLine!="[Data]");
        getline(infs, aLine);   // the \n after [Data] line
        getline(infs, aLine);   // the header line of [Data] section

        // going through the data section
        string prevID ("--"), curAnimID;
        unsigned curIDindex(0);
        string lineGenotype;
        lineGenotype.reserve(numberSNPs*3);     // reserve space for efficiency

        // write out the new block
//        for(i=0; i<vecReservedBlock.size(); ++i)
//            opfs << vecReservedBlock[i] << endl;
        // write out the full pathname of the input finalReport file
        opfs << "#" << finalReportFullPath;
        for(i=0; i<vecSampleID.size(); ++i)
            opfs << "#" << vecSampleID[i];
        opfs << endl;
        progClear();

        while(getline(infs, dataLine))
        {
            sx=dataLine.split('\t');
            if(prevID!=sx[1]) // new Animal
            {
                prevID=sx[1];
                ++curIDindex;
                progShow('W', 100*curIDindex/numberSamples);

                if(curIDindex>1)// && vecSampleID2bAppended[ind]==curAnimID)// write
                {
//                    ++ind;
                    opfs
                    << left << setw(32) << curAnimID
                    << lineGenotype << endl;
                }
                // renew vars
                lineGenotype="";
                curAnimID=sx[1];
                // replace whitespaces with '_'
                if(curAnimID.find(" ")!=string::npos)
                    for(i=0; i<curAnimID.length(); ++i)
                        if(curAnimID[i]==' ') curAnimID[i]='_';
                if(mapXref.find(curAnimID)!=mapXref.end())
                    curAnimID=mapXref[curAnimID];
            }
            lineGenotype+=sx[posAllele1AB];
            lineGenotype+=sx[posAllele2AB];
            lineGenotype+=" ";
        }
        // write out the data of the last animal
//        if(vecSampleID2bAppended[ind]==curAnimID)
            opfs
            << left << setw(32) << curAnimID
            << lineGenotype << endl;

        infs.close();
        opfs.close();

        /* rename */

        i=rename(tempFile.c_str(), ABgenotypeFileName.c_str());
        progClear();
        if(i==0)
        {
            I("sisssis", 2, "Finally,", vecSampleID.size(), "new sample(s) BLOCK-REPLACED to [", ABgenotypeFileName.c_str(), "] within",
              (clock()-startTime)/CLOCKS_PER_SEC, "seconds. √");
            return true;
        }
        // else
        E("ss", "ERROR: Could not rename the temp file to the AB-genotyped file:", tempFile.c_str());
        return false;
    }

    void finalReport::clear()
    {
        vector<string>().swap(vecSampleID);
        vector<string>().swap(vecSNPName);
        vector<string>().swap(vecStrx);
//        vector<string>().swap(vecSampleID2bAppended);
//        set<string>().swap(setABSampleID);
        map<string, string>().swap(mapXref);
    }
}
