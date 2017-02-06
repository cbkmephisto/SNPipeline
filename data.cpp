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
//  data.cpp
//  SNPipeline
//
//  Created by Hailin SU on 1/29/14.
//

#include "data.h"
#include <cstdio>
#include <climits>
#include <iomanip>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)

using namespace std;

namespace libcbk
{
    /***************************************************************************pipeline
     *
     *  runs the whole routine
     */
    void Data::pipeline(const string fileName)
    {
        startTime = clock();
        clear();
        static int ret;
        ret=preCheck(fileName);
        //  > 0, go ahead
        //  < 0, ERROR
        if(ret<0)
            return;

        ret=postCheck();
        //  == 0, info already exists, continue
        //   < 0, ERROR
        if(ret<0)
            return;
        else if(ret==1) // simple append
            DT2AB_simpleAppend();
        else            // merging append
            DT2AB_mergingAppend();
        clear();
    }

    /***************************************************************************checkAffIll
     *
     *  checks if the format of given file is
     *      affymetrix input
     *  or
     *      illumina finalReport
     *  stores in class member strAffIll
     */
    void Data::checkAffIll(const string fileName)
    {
        ifstream infsx (fileName.c_str());
        if(!infsx)
        {
            E("ss", "File not found:", fileName.c_str());
            exit(-1001);
        }
        string aLine;
        infsx >> aLine;  // header: [Header]==Illumina, probeset_id==Affymetrix
        infsx.close();
        if(aLine=="[Header]")
            strAffIll="Ill";
        else if(aLine=="probeset_id")
        {
            strAffIll="Aff";
            content="Affymetrix";
        }
        else
        {
            E("ss", "File is neither an Illumina finalReport file nor an Affymetrix genotype file:", fileName.c_str());
            exit(-1002);
        }
    }

    /***************************************************************************clear
     *
     *  clean working space
     */
    void Data::clear()
    {
        strAffIll="";
        vector<string>().swap(vecSampleID);
        vector<string>().swap(vecSNPName);
        vector<string>().swap(vecStrx);
        vector<unsigned char>().swap(vecGenotypes);

        map<string, string>().swap(mapXref);
    }

    /***************************************************************************preCheck for all
     *
     *  clean working space
     */
    int Data::preCheck(const string fileName)
    {
        // ******** init static map/vector
        /*
        if(mapAB2uc.find("AB")==mapAB2uc.end())
        {
            mapAB2uc["NoCall"]=0;
            mapAB2uc["--"]=0;
            mapAB2uc["AA"]=1;
            mapAB2uc["BB"]=2;
            mapAB2uc["AB"]=3;
            mapAB2uc["BA"]=3;
        }
        // * newer algorithm for coding AB before above
         * 'A'=65, 'B'=66, 'N'=78, 'o'=111, '-'=45
         * "AA"='A'+'A'=130
         * "BB"='B'+'B'=132
         * "AB"='A'+'B'=131
         * "No"='N'+'o'=189
         * "--"='-'+'-'=90
         */
        if(mapAB2uc.find(130)==mapAB2uc.end())// at last I didn't use this map, just leave it here
        {
            mapAB2uc[189]=mapAB2uc[90]=0;   // NoCall --
            mapAB2uc[130]=1;    // AA
            mapAB2uc[132]=2;    // BB
            mapAB2uc[131]=3;    // AB
        }
        if(vecUC2STR.size()==0)
        {
            vecUC2STR.push_back("--");// 0
            vecUC2STR.push_back("AA");// 1
            vecUC2STR.push_back("BB");// 2
            vecUC2STR.push_back("AB");// 3
        }

        // ******** check input file type
        checkAffIll(fileName);

        // ******** load the data
        if(strAffIll=="Aff")
        {
            if(loadA(fileName)<0)
                return -1;
        }
        else
            if(loadI(fileName)<0)
                return -1;

        unsigned i, j;

        // ******** replace whitespaces with '_'
        for(i=0; i<numberSNPs; ++i)
            if(vecSNPName[i].find(" ")!=string::npos || vecSNPName[i].find("#")!=string::npos)
            {
                for(j=0; j<vecSNPName[i].length(); ++j)
                    if(vecSNPName[i][j]==' ' || vecSNPName[i][j]=='#')
                        vecSNPName[i][j]='_';
            }

        for(i=0; i<numberSamples; ++i)
            if(vecSampleID[i].find(" ")!=string::npos || vecSampleID[i].find("#")!=string::npos)
            {
                for(j=0; j<vecSampleID[i].length(); ++j)
                    if(vecSampleID[i][j]==' ' || vecSampleID[i][j]=='#')
                        vecSampleID[i][j]='_';
            }

        // ******** storing extra info into the object

//        char buff[512];
//        realpath(fileName.c_str(), buff);
        char* buff=realpath(fileName.c_str(), NULL);
        for(i=0; i<512; ++i)
        {
            if(buff[i]=='/')
                j=i;    // record the last /
//            else if(buff[i]==' ')
//                buff[i]='_';
            else if(buff[i]==0)
            {
                buff[j+1] = 0;  // discard fileName
                break;
            }
        }
        fullPath=buff;

        string stmp;
        stringstream sst;
        sst << numberSNPs;
        sst >> stmp;
        sst.clear();
        ABgenotypeFileName      = "ab-genotyped-" + stmp + "-" + content;
        ABgenotypeXrefFileName  = fullPath + "xref";
        free(buff);
        return 1;
    }

    /***************************************************************************loadA for affy
     *
     * read basic info from affyInput file, such as
     *      numberSNPs
     *      numberSamples
     *      vecSNPIDs
     *      vecSampleIDs
     *      vecGenotypes
     * and checks for duplicates
     */
    int Data::loadA(const string fileName)
    {
        //        I("s", 2, "Pre-checking...");
        infs.open(fileName.c_str());
        if(!infs)
        {
            E("ss", "File not found:", fileName.c_str());
            return -1;
        }

        inputFileName=fileName;
        numberSNPs = numberSamples = 0;
        string          stmp;
        splitablestring aLine;
        stringstream    sst;
        //        vector<string>  vecStrx;
        unsigned        i;

        const unsigned long long FILE_SIZE (0.01*tellFileSize(fileName));
        // ****************** clear working space
        vecSNPName.reserve(650000); // reserve 650k
        // ****************** read the whole file
        // 1st line: the header
        if(getline(infs, aLine))
        {
            vecStrx = aLine.split('\t');
            // vecStrx[0]==the 'probeset_id' flag
            if(vecStrx.size()/2==(1+vecStrx.size())/2)
            {
                E("s", "Wrong file format: column number != 1 + 2 * numberSampleIDs. Aborting.");
                infs.close();
                exit(-999);
            }

            for(i=1; i<vecStrx.size(); ++i)
            {
                //                sst >> stmp;   // each sample ID appears twice
                if(vecStrx[i]!=vecStrx[i+1])
                {
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
                    E("ss", "Aborted due to duplicated Sample ID found:", stmp.c_str());
                    infs.close();
                    exit(-1006);
                }
            }

            set<string>().swap(setSampleIDs);
            numberSamples=vecSampleID.size();
            for(i=0; i<vecSampleID.size(); ++i)
            {
                for(unsigned j=0; j<vecSampleID[i].size(); ++j)
                    if(vecSampleID[i][j]=='#')
                        vecSampleID[i][j]='_';
            }
        }
        //        while (infs>>stmp)

        unsigned curLine(1), curCharCount(0);
        unsigned char ctmp(0), cgtmp;
        vecGenotypes.reserve(650000*(1+numberSamples)/4);   // reserve space

        // rest of the file: genotype data
        while(getline(infs, aLine))
        {
            if(++curLine%6144==2) // about 1% of the 650k SNPs
                progShow('L', unsigned(infs.tellg()/FILE_SIZE));
            vecStrx = aLine.split('\t');
            if(vecStrx.size()%2==0)    // could not get the 2nd column
            {
                cout << endl;
                E("sssi", "Aborted due to bad data line (could not get the 2nd column for a sample) for SNP:", vecSNPName.end()->c_str(),
                  "in line", curLine);
                infs.close();
                exit(-1008);
            }
            vecSNPName.push_back(vecStrx[0]);

            for(i=1; i<vecStrx.size(); i+=2)   // vecStrx[i] is the AB genotype
            {
                cgtmp=vecStrx[i][0]+vecStrx[i][1];
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
                /* newer algorithm for coding AB before above
                 * 'A'=65, 'B'=66, 'N'=78, 'o'=111, '-'=45
                 * "AA"='A'+'A'=130
                 * "BB"='B'+'B'=132
                 * "AB"='A'+'B'=131
                 * "No"='N'+'o'=189
                 * "--"='-'+'-'=90
                 mapAB2uc[189]=mapAB2uc[90]=0;   // NoCall --
                 mapAB2uc[130]=1;    // AA
                 mapAB2uc[132]=2;    // BB
                 mapAB2uc[131]=3;    // AB
                if(mapAB2uc.find(cgtmp)==mapAB2uc.end())    // not a valid genotype
                {
                    cout << endl;
                    E("sssi", "Aborted due to invalid SNP genotype:", vecStrx[i].c_str(), "in line", curLine);
                    infs.close();
                    exit(-1018);
                }
                */
                ctmp<<=2;
                switch(cgtmp)
                {
                    case 131:   // AB
                        cgtmp=3;
                        break;
                    case 130:   // AA
                        cgtmp=1;
                        break;
                    case 132:   // BB
                        cgtmp=2;
                        break;
//                    case 90:    // --
                    case 189:   // NoCall
                        cgtmp=0;
                        break;
                    default:
                        cout << endl;
                        E("sssi", "Aborted due to invalid SNP genotype:", vecStrx[i].c_str(), "in line", curLine);
                        infs.close();
                        exit(-1018);
                }
                //ctmp+=mapAB2uc[cgtmp];   // left shift 2 bits and store current genotype at the lower 2 bits
                ctmp+=cgtmp;
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
        infs.close();
        progClear();
        numberSNPs=vecSNPName.size();

        // ******** checks for duplicated SNPNames
        set<string>().swap(setSNPIDs);
        for(i=0; i<vecSNPName.size(); ++i)
            setSNPIDs.insert(vecSNPName[i]);
        if(vecSNPName.size()!=setSNPIDs.size())
        {
            E("s", "Duplicated SNP name(s) found:");
            unsigned j;
            for(i=0; i<vecSNPName.size()-1; ++i)
                for(j=i+1; j<vecSNPName.size(); ++j)
                {
                    if(vecSNPName[i]==vecSNPName[j])
                    {
                        cout << "\t\t\t" << vecSNPName[i] << endl;
                        break;  // break the inner loop
                    }
                }
            E("s", "Aborting... X");
            return -8;
        }
        set<string>().swap(setSNPIDs);
        return 1;
    }

    /***************************************************************************loadI for Illu
     *
     * read basic info from affyInput file, such as
     *      numberSNPs
     *      numberSamples
     *      content
     *      vecSNPIDs
     *      vecSampleIDs
     *      vecGenotypes
     *      posSNPName
     *      posSampleID
     *      posAllele1AB
     *      posAllele2AB
     * and checks for duplicates
     */
    int Data::loadI(const string fileName)
    {
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
        unsigned i(0), j;
        string aLine;

        // **************** [Header] part
        while(infs >> aLine && aLine!="[Data]")
        {
            //            if(aLine=="[Header]")
            //                I("s", 2, aLine.c_str());
            if(++i>200) break;       // avoid large file
            if(aLine=="Content")
            {
                //infs >> content;
                getline(infs, content);
                for(unsigned c=0; c<content.length(); ++c)
                {
                    if(content[c]==' ' || content[c]==9) // 9=ASCII(TAB)
                       content[c]='_';
                }
                // 20150901
                while(content[0]=='_')
                {
                    if(content.length()>1 )    // prevent content starting with underscore '_'
                        content=content.substr(1);
                }
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
        unsigned noS1p (0.01*numberSamples);
        if(noS1p==0)    // avoid divided by zero
            noS1p=1;
        vecGenotypes.reserve(numberSamples*numberSNPs/4);
        unsigned char ctmp(0), cgtmp;
        unsigned curCharCount(0);

        // we choose to believe H0: numberSamples and numberSNPs from the header part

        //// ****************  read 1st sample
        progShow('L', 0);
        for(j=0; j<numberSNPs; ++j)
        {
            if(getline(infs, dataLine))
            {
                vecStrx=dataLine.split('\t');
                setSNPIDs.insert(vecStrx[0]);
                vecSNPName.push_back(vecStrx[0]);
                setSampleIDs.insert(vecStrx[1]);
                // encoding, detailed see loadA
                cgtmp=vecStrx[posAllele1AB][0]+vecStrx[posAllele2AB][0];
                ctmp<<=2;
                switch(cgtmp)
                {
                    case 131:   // AB
                        cgtmp=3;
                        break;
                    case 130:   // AA
                        cgtmp=1;
                        break;
                    case 132:   // BB
                        cgtmp=2;
                        break;
                    case 90:    // --
                        //                case 189:   // NoCall
                        cgtmp=0;
                        break;
                    default:
                        cout << endl;
                        E("sssss", "Aborted due to invalid SNP genotype:", (vecStrx[posAllele1AB]+vecStrx[posAllele2AB]).c_str(), "in line", vecStrx[0].c_str(), vecStrx[1].c_str());
                        infs.close();
                        exit(-1018);
                }
                ctmp+=cgtmp;
                if(++curCharCount==4) // full, plug in and refresh
                {
                    vecGenotypes.push_back(ctmp);
                    ctmp=curCharCount=0;
                }
            }
            else
            {
                cout << endl;
                E("sisis", "EOF while loading 1st sample, only", j, "SNPs read but numberSNPs =", numberSNPs, "from [Header] part.");
                E("s", "BBBBBBBBBBBBBBBBBBBad finalReport file, aborting. X");
                infs.close();
                return -6;
            }
        }

        //// ****************  check 1st sample

        // checks for duplicated SNP Names
        if(vecSNPName.size()!=setSNPIDs.size())
        {
            cout << endl;
            E("s", "Duplicated SNP name(s) found:");
            for(i=0; i<vecSNPName.size()-1; ++i)
                for(j=i+1; j<vecSNPName.size(); ++j)
                {
                    if(vecSNPName[i]==vecSNPName[j])
                    {
                        cout << "\t\t\t" << vecSNPName[i] << endl;
                        break;  // break the inner loop
                    }
                }
            E("s", "Aborting... X");
            return -8;
        }
        set<string>().swap(setSNPIDs);
        if(vecSNPName.size()!=numberSNPs)
        {
            I("sis", 3, "number of SNPNames read = ", vecSNPName.size(), "X");
            E("sis", "ERROR: SNP_Names_read != numberSNPs", numberSNPs, "from [Header] section, bad finalReport file, aborting. X");
            return -20;
        }

        // checks for 1st Sample ID
        if(setSampleIDs.size()!=1)
        {
            cout << endl;
            E("sis", "ERROR while loading 1st sample,", setSampleIDs.size(), "sampleIDs read in the 1st sample block according to the header info.");
            E("s", "BBBBBBBBBBBBBBBBBBBad finalReport file, aborting. X");
            infs.close();
            return -9;
        }
        vecSampleID.push_back(*setSampleIDs.begin());

        //// ****************  read other samples
        vector<string> vecTmpSNPNames;
        set<string> setTmpSampleIDs;

        for(i=1; i<numberSamples; ++i)
        {
            if(i%noS1p==0)
            {
                if(noS1p==1)
                    progShow('L', i*100/numberSamples);
                else
                    progShow('L', i/noS1p);
            }
            for(j=0; j<numberSNPs; ++j)
            {
                if(getline(infs, dataLine))
                {
                    vecStrx=dataLine.split('\t');
                    vecTmpSNPNames.push_back(vecStrx[0]);
                    setTmpSampleIDs.insert(vecStrx[1]);
                    // encoding, detailed see loadA
                    cgtmp=vecStrx[posAllele1AB][0]+vecStrx[posAllele2AB][0];
                    ctmp<<=2;
                    switch(cgtmp)
                    {
                        case 131:   // AB
                            cgtmp=3;
                            break;
                        case 130:   // AA
                            cgtmp=1;
                            break;
                        case 132:   // BB
                            cgtmp=2;
                            break;
                        case 90:    // --
                            //                case 189:   // NoCall
                            cgtmp=0;
                            break;
                        default:
                            cout << endl;
                            E("sssss", "Aborted due to invalid SNP genotype:", (vecStrx[posAllele1AB]+vecStrx[posAllele2AB]).c_str(), "in line", vecStrx[0].c_str(), vecStrx[1].c_str());
                            infs.close();
                            exit(-1018);
                    }
                    ctmp+=cgtmp;
                    if(++curCharCount==4) // full, plug in and refresh
                    {
                        vecGenotypes.push_back(ctmp);
                        ctmp=curCharCount=0;
                    }
                }   // end of if not EOF
                else
                {
                    cout << endl;
                    E("siss", "UNEXPECTED EOF while loading data for sample", i+1, ":", vecStrx[1].c_str());
                    E("si", "    sample offset (0-based):", i);
                    E("si", "    SNP_ID offset (0-based):", j);
                    E("s", "One possible cause of this error, is that in the final report file, the number of samples reported in the header part, is larger than the number of samples it is actually holding. Try launching the following unix command to see actual number of samples in the final report file:");
                    E("s", "    awk 'NR>10{print $2}' <finalReportFile> | uniq -c | awk 'NF==2' | wc -l");
                    E("s", "BBBBBBBBBBBBBBBBBBBad finalReport file, aborting. X");
                    infs.close();
                    return -6;
                }   // end of if not EOF
            }   // end of for each SNP

            //// ****************  check ith sample
            // checks for ist Sample ID
            if(setTmpSampleIDs.size()!=1)
            {
                cout << endl;
                E("sisis", "ERROR:", setTmpSampleIDs.size(), "sampleIDs read in the", i+1, "th sample block according to the header info:");
                for(set<string>::iterator sit = setTmpSampleIDs.begin(); sit!=setTmpSampleIDs.end(); ++sit)
                    cerr << "\t" << (*sit);
                cerr << endl;
                E("s", "BBBBBBBBBBBBBBBBBBBad finalReport file, aborting. X");
                infs.close();
                return -19;
            }
            vecSampleID.push_back(*setTmpSampleIDs.begin());
            setSampleIDs.insert(*setTmpSampleIDs.begin());
            set<string>().swap(setTmpSampleIDs);

            // checks for orders of SNP Names
            for(j=0; j<numberSNPs; ++j)
                if(vecTmpSNPNames[j]!=vecSNPName[j])
                {
                    cout << endl;
                    infs.close();
                    E("sssssss", "ERROR: different SNP names/order found, Sample [", vecSampleID[vecSampleID.size()-1].c_str(),
                      "] - SNP [", vecTmpSNPNames[j].c_str(), "], SNPName should be {", vecSNPName[j].c_str(), "}, bad finalReport file, aborting. X");
                    infs.close();
                    return -8;
                }
            vector<string>().swap(vecTmpSNPNames);
        }   // end of for each sample

        if(curCharCount!=0) // somthing-left non-full
        {
            // move them to the left
            ctmp <<= 2*(4-curCharCount);
            vecGenotypes.push_back(ctmp);
        }

        infs.close();
        progClear();

        if(vecSampleID.size()!=numberSamples)
        {
            I("sis", 3, "number of SampleIDs read = ", vecSampleID.size(), "X");
            E("sis", "ERROR: Sample_IDs_read != numberSamples", numberSamples, "from [Header] section, bad finalReport file, aborting. X");
        }
        getline(infs, dataLine);// test
        if(!infs.eof() || dataLine.length()!=0)
        {
            E("s", "ERROR: I thought I've read all the data (numberSamples*numberSNPs entries), but there's still something left in the file.");
            E("s", "Please re-check if numberSamples and numberSNPs in the header part is consistent with the data body. Bad finalReport file, aborting. X");
            cerr << dataLine << endl;
            exit (-9);
        }

        if(vecSampleID.size()!=setSampleIDs.size())
        {
            I("s", 1, "Duplicated sampleID(s) found:");
            for(i=0; i<vecSampleID.size()-1; ++i)
                for(j=i+1; j<vecSampleID.size(); ++j)
                {
                    if(vecSampleID[i]==vecSampleID[j])
                    {
                        cout << "\t\t\t" << vecSampleID[i] << endl;
                        break;  // break the inner loop
                    }
                }
//            E("s", "Aborting... X");
//            return -28;
            I("s", 1, "Try continuing with duplicated sampleIDs"); // 20150812
        }
        set<string>().swap(setSampleIDs);

        return 1;
    }


    /***************************************************************************postCheck
     *
     * basic info related to AB-genotype file
     */
    int Data::postCheck()
    {
        //        I("s", 2, "Post-checking...");
        splitablestring     aLine;
        stringstream        ss;
        vector<string>      sx;
        string              stmp;
        bool                bXref(true);
        unsigned            i(0);

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
//                if(sx.size()!=2)        // 20150813: allowing more than 2 cols, xref [0] -> [last]
                if(sx.size()<2)
                {
                    cout << endl;
                    // wrong xref file, skipping
                    E("s", "WARNING: illegal xref file read (not a >=2-column-whole-through file), skipping xref functions. #");
                    map<string, string>().swap(mapXref);
                    bXref=false;
                    break;
                }
//                mapXref[sx[0]]=sx[1];
                mapXref[sx[0]]=sx[sx.size()-1];
                vector<string>().swap(sx);
            }
            infs.close();

            if(bXref)   // update the vecSampleID
            {
                for(i=0; i<vecSampleID.size(); ++i)
                    if(mapXref.find(vecSampleID[i])!=mapXref.end()) // got a xref
                    {
                        vecSampleID[i]=mapXref[vecSampleID[i]];
                    }
            }
        }

        // **************** checks for duplication in both sample IDs and SNP IDs
        for(i=0; i<vecSampleID.size(); ++i)
            setSampleIDs.insert(vecSampleID[i]);

        if(vecSampleID.size()!=setSampleIDs.size())
        {
            cout << endl;
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

        infs.open(ABgenotypeFileName.c_str());
        if(!infs)
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
                cout << endl;
                infs.close();
                E("s", "The existing AB-genotype file has a wrong file format. Aborting. X");
                return -1;
            }
            // **************** todo: continue with the SNP_IDs
            progShow('P', 75+clock()%11);
            for(i=1; i<sx.size(); ++i)
            {
                if(sx[i]!=vecSNPName[i-1])
                {
                    cout << endl;
                    E("sssss", "ERROR: different SNP names/order found, SNP", sx[i].c_str(), "in existed ab-genotype file, while", vecSNPName[i-1].c_str(), "in the same position of input file, aborting. X");
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

    /***************************************************************************write
     *
     *  write: distribute to writeA/writeI according to strAffIll
     */
    void Data::write(ofstream &ofs)
    {
        static unsigned i;
        static string stmp;
        //        static unsigned char ctmp;

        if(!abFileExists)
        {
            // write out the SNP names for the whole file
            ofs << "#SNP_IDs";
            for(i=0; i<numberSNPs; ++i)
                ofs << "#" << vecSNPName[i];
            ofs << endl;
            abFileExists=true;
        }
        // write out the full pathname of the input file
        ofs << "#" << fullPath;
        for(i=0; i<numberSamples; ++i)
            ofs << "#" << vecSampleID[i];
        ofs << endl;
        
        // write out the genotype
        if(strAffIll=="Aff")
            writeA(ofs);
        else
            writeI(ofs);
    }

    /***************************************************************************writeA
     *
     * write vecGenotypes to given ofstream: Affy
     */
    void Data::writeA(ofstream &ofs)
    {
        static unsigned i, j, shift;
        unsigned long segIndex(0);       // each char for 4 segments [ 11 01 00 10 ]
        //        unsigned long charIndex(0);      // the index of vecGenotype for the current Animal & current SNP
        //        unsigned char mask(0);      // the mask for current Animal for the current vecGenotypes char
        const unsigned long numx(numberSNPs*numberSamples);
        for(i=0; i<numberSamples; ++i)
        {
            if(i%7==0)  // show % every 7 samples
                progShow('W', 100*i/numberSamples);

            ofs << left << setw(32) << vecSampleID[i] << " ";   // 20150812: prevent too long ID
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
                shift = 2 * (3 - (segIndex % 4));
                ofs << vecUC2STR[( (3 << shift) & vecGenotypes[segIndex/4] ) >> shift] << " ";
            }
            ofs << endl;
        }
        progClear();
    }

    /***************************************************************************writeI
     *
     * write vecGenotypes to given ofstream: Illu
     */
    void Data::writeI(ofstream &ofs)
    {
        static unsigned i, j;
        unsigned charIndex (0);
        unsigned long itmp;
        for(i=0; i<numberSamples; ++i)
        {
            if(i%7==0)  // show % every 7 samples
                progShow('W', 100*i/numberSamples);

            ofs << left << setw(32) << vecSampleID[i] << " ";   // 20150812: prevent too long ID
            itmp=i*numberSNPs;
            for(j=0; j<numberSNPs; ++j)
            {
                switch((itmp+j)%4)
                {
                    case 0:
                        ofs << vecUC2STR[(vecGenotypes[charIndex] & 192) >> 6] << " "; // 0b11000000
                        break;
                    case 1:
                        ofs << vecUC2STR[(vecGenotypes[charIndex] & 48) >> 4] << " "; // 0b00110000
                        break;
                    case 2:
                        ofs << vecUC2STR[(vecGenotypes[charIndex] & 12) >> 2] << " "; //0b00001100
                        break;
                    case 3:
                        ofs << vecUC2STR[vecGenotypes[charIndex++] & 3] << " "; // 0b00000011
                }
            }
            ofs << endl;
        }
        progClear();
    }

    /***************************************************************************simpleAppend
     *
     * input file proceed to AB-genotype file with simply appending the content
     */
    bool Data::DT2AB_simpleAppend()
    {
        // ****************** output
        ofstream ofs(ABgenotypeFileName.c_str(), ofstream::app | ofstream::out);  // appending mode
        if(!ofs)
        {
            cout << endl;
            E("s", "Could not write to file:", ABgenotypeFileName.c_str());
            return false;
        }

        write(ofs);

        ofs.close();
        I("sisssis", 2, "Finally,", numberSamples, "new sample(s) WRITEN into [", ABgenotypeFileName.c_str(), "] within",
          (clock()-startTime)/CLOCKS_PER_SEC, "seconds. √");
        return true;
    }

    /***************************************************************************mergingAppend
     *
     * input file proceed to AB-genotype file with merging appending the content
     */
    bool Data::DT2AB_mergingAppend()
    {
        // creating a temporary ab-genotype file without the shared block
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

}