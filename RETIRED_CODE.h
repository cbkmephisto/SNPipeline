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
//  RETIRED_CODE.h
//  SNPipeline
//
//  Created by Hailin SU on 1/28/14.
//

#ifndef SNPipeline_RETIRED_CODE_h
#define SNPipeline_RETIRED_CODE_h

//**************************************************************************************************************************
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
int finalReport::preCheckA(const string fileName)
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
    unsigned i (0);
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
    string curAnimID, prevID("--"), sidtmp;
    unsigned curIDindex(0), SNPIDx(0);
    while(getline(infs, dataLine))
    {
        vecStrx=dataLine.split('\t');
        if(prevID!=vecStrx[1]) // new Animal posSampleID==1
        {
            prevID=vecStrx[1];
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
                    //                        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << flush;
                    //                        I("si", 3, "SNP IDs read  = ", vecSNPName.size());
                }
            }
            // push sample id
            curAnimID=vecStrx[1];
            // replace whitespaces with '_'
            if(curAnimID.find(" ")!=string::npos)
                for(i=0; i<curAnimID.length(); ++i)
                    if(curAnimID[i]==' ') curAnimID[i]='_';
            vecSampleID.push_back(curAnimID);
        }
        if(curIDindex==1)   // data of 1st sample: fill SNP Names
        {
            // replace whitespaces with '_'
            sidtmp=vecStrx[0];
            if(sidtmp.find(" ")!=string::npos)
                for(i=0; i<sidtmp.length(); ++i)
                    if(sidtmp[i]==' ') sidtmp[i]='_';
            vecSNPName.push_back(sidtmp);
        }
        else                // data of others: check the SNP IDs
        {
            if(vecSNPName[SNPIDx++]==vecStrx[0])
                continue;
            progClear();
            E("sssssss", "ERROR: different SNP names/order found, Sample", curAnimID.c_str(),
              "- SNP", vecStrx[posSNPName].c_str(), "should be", vecSNPName[SNPIDx-1].c_str(), ", bad finalReport file, aborting. X");
            return -6;
        }
    }

    progClear();
    infs.close();
    if(vecSampleID.size()!=numberSamples)
    {
        I("sis", 3, "SampleID read = ", vecSampleID.size(), "X");
        E("s", "ERROR: Sample_IDs_read != numberSamples from [Header] section, bad finalReport file, aborting. X");
        return -7;
    }
    //        I("si", 3, "SampleID read = ", vecSampleID.size());

    //        I("s", 3, "Checking for duplicated IDs...");

    // **************** checks for duplication in both sample IDs and SNP IDs
    set<string> setSampleIDs, setSNPNames;
    for(i=0; i<vecSampleID.size(); ++i)
        setSampleIDs.insert(vecSampleID[i]);
    if(vecSampleID.size()!=setSampleIDs.size())
    {
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
    unsigned j(0);
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
    return 1;
}

//**************************************************************************************************************************
// load affyInput file into vecGenotypes
//  DEPRECATED: IMPLEMENTED IN PRECHECK() Jan 27 2014
void affyInput::load()
{
    // ****************** read the whole file
    if(mapAB2uc.find("AB")==mapAB2uc.end()) // need init
    {
        mapAB2uc["NoCall"]=0;
        mapAB2uc["--"]=0;
        mapAB2uc["AA"]=1;
        mapAB2uc["BB"]=2;
        mapAB2uc["AB"]=3;
        mapAB2uc["BA"]=3;
    }
    infs.open(inputFileName.c_str());
    string aLine, stmp;
    stringstream sst;
    // 1st line: the header
    getline(infs, aLine);
    unsigned curp(0), lastp(0), curLine(1), curCharCount(0);
    unsigned char ctmp(0);
    const unsigned long long FILE_SIZE (tellFileSize(inputFileName));
    vecGenotypes.reserve(numberSNPs*(1+numberSamples)/4);   // reserve space

    // rest of the file: data
    while (getline(infs, aLine))
    {
        ++curLine;
        sst << aLine;
        if(sst >> stmp)   // 1st column is the SNP ID
        {
            // the genotypes
            while(sst>>aLine)   // aLine is AB genotype
            {
                if(!(sst>>stmp))    // do not get the 2nd column
                {
                    progClear();
                    E("sssi", "Aborted due to bad data line (could not get the 2nd column for a sample) for SNP:", vecSNPName.end()->c_str(),
                      "in line", curLine);
                    infs.close();
                    exit(-1008);
                }
                if(mapAB2uc.find(aLine)==mapAB2uc.end())    // not a valid genotype
                {
                    progClear();
                    E("sssi", "Aborted due to invalid SNP genotype:", aLine.c_str(), "in line", curLine);
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
                ctmp<<=2;   // left shift 2 bits
                ctmp+=mapAB2uc[aLine];
                curCharCount++;
                if(curCharCount==4) // full, plug in and refresh
                {
                    vecGenotypes.push_back(ctmp);
                    ctmp=curCharCount=0;
                }
            }
        }
        sst.clear();
        curp=100.0*infs.tellg()/FILE_SIZE;
        if(curp!=lastp)
        {
            lastp=curp;
            progShow('L', lastp);
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
}

//**************************************************************************************************************************
// memory consuming algorithm
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
    unsigned i (0);
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
    // reads all the raw-Sample/SNP IDs
    vector<string> vecRawSampleIDs;
    vector<string> vecRawSNPIDs;
    const unsigned long long FILE_SIZE (tellFileSize(fileName));
    const unsigned long vecSize(numberSamples*numberSNPs);
    unsigned p1p     (0.015*vecSize);
    vecRawSampleIDs.reserve(vecSize);
    vecRawSNPIDs.reserve(vecSize);
    unsigned lastp(0);
    unsigned char ctmp(0);
    string stmp;
    unsigned curCharCount(0);

    if(p1p<2)
        p1p=2;
    while(getline(infs, dataLine))
    {
        if(++lastp%p1p==0)
            progShow('L', (unsigned)(95*infs.tellg()/FILE_SIZE));   // max 95%
        vecStrx=dataLine.split('\t');
        vecRawSNPIDs.push_back(vecStrx[0]);
        vecRawSampleIDs.push_back(vecStrx[1]);
        // encoding
        stmp=vecStrx[posAllele1AB]+vecStrx[posAllele2AB];
        if(mapAB2uc.find(stmp)==mapAB2uc.end())    // not a valid genotype
        {
            cout << endl;
            E("sssss", "Aborted due to invalid SNP genotype:", stmp.c_str(), "in line", vecStrx[0].c_str(), vecStrx[1].c_str());
            infs.close();
            exit(-1018);
        }
        /* new algorithm for storing genotype vector in the memory
         *
         *  AA    BB    AB/BA NoCall/--
         *  1     2     3     0
         *
         *  each unsigned char holds up to 4 SNPs
         *
         *  unsigned char (5) = 0000 0101 = -- -- AA AA
         *  unsigned char (3) = 0000 0011 = -- -- -- AB
         *  unsigned char (?) = 1110 0100 = AB BB AA --
         */
        ctmp = (ctmp << 2) + mapAB2uc[stmp];   // left shift 2 bits and store current genotype at the lower 2 bits
        //                ctmp+=mapAB2uc[aLine];
        if(++curCharCount==4) // full, plug in and refresh
        {
            vecGenotypes.push_back(ctmp);
            ctmp=curCharCount=0;
        }
    }

    if(curCharCount!=0) // somthing-left non-full
    {
        // move them to the left
        ctmp <<= 2*(4-curCharCount);
        vecGenotypes.push_back(ctmp);
    }

    infs.close();

    if(vecRawSNPIDs.size()!=vecSize)    // check if nSNPs and nSamples consistent with file data body
    {
        cout << endl;
        E("s", "ERROR: potentially index out of boundary, please re-check if numberSamples and numberSNPs in the header part is consistent with the data body. Bad finalReport file, aborting. X");
        I("si", 1, "[Header] nSNP * nSample =", vecSize);
        I("si", 1, "but tot data entry read =", vecRawSampleIDs.size());
        return -15;
    }
    progClear();

    // filling vecSNPName and vecSampleID
    string prevID("--"), sidtmp;
    unsigned curIDindex(0), SNPIDx(0);
    for(i=0; i<vecSize; ++i)
    {
        if(prevID!=vecRawSampleIDs[i]) // new Animal
        {
            prevID=vecRawSampleIDs[i];

            progShow('L', 95+5*curIDindex/numberSamples);
            ++curIDindex;

            if(curIDindex>1)// when SNP name reading has been done,
            {
                SNPIDx=0;
                if(curIDindex==2)   // collect the info
                {
                    if(vecSNPName.size()!=numberSNPs)
                    {
                        cout << endl;
                        I("sis", 3, "SNP IDs read  = ", vecSNPName.size(), "X");
                        E("sis", "ERROR: SNP_IDs_read != numberSNPs (", numberSNPs, ") from [Header] section, bad finalReport file, aborting. X");
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
            if(SNPIDx<vecSNPName.size() && vecSNPName[SNPIDx++]==vecRawSNPIDs[i])
                continue;
            cout << endl;
            if(SNPIDx==vecSNPName.size())
                E("s", "ERROR: index out of boundary for nSNPs per animal, please re-check the inputed bad finalReport file, aborting. X");
            else
                E("sssssss", "ERROR: different SNP names/order found, Sample [", prevID.c_str(),
                  "] - SNP [", vecRawSNPIDs[i].c_str(), "], SNPName should be {", vecSNPName[SNPIDx-1].c_str(), "}, bad finalReport file, aborting. X");
            return -6;
        }
    }

    vector<string>().swap(vecRawSNPIDs);
    vector<string>().swap(vecRawSampleIDs);

    if(vecSampleID.size()!=numberSamples)
    {
        cout << endl;
        I("sis", 3, "number of SampleIDs read = ", vecSampleID.size(), "X");
        E("sis", "ERROR: Sample_IDs_read != numberSamples", numberSamples, "from [Header] section, bad finalReport file, aborting. X");
        return -7;
    }
    if(vecSNPName.size()!=numberSNPs)
    {
        cout << endl;
        I("sis", 3, "number of SNPNames read = ", vecSNPName.size(), "X");
        E("sis", "ERROR: SNP_Names_read != numberSNPs", numberSNPs, "from [Header] section, bad finalReport file, aborting. X");
        return -7;
    }

    for(i=0; i<vecSampleID.size(); ++i)
    {
        for(unsigned j=0; j<vecSampleID[i].size(); ++j)
            if(vecSampleID[i][j]=='#')
                vecSampleID[i][j]='_';
    }
    progClear();
    return 1;
}
//**************************************************************************************************************************
// slow checking Sample/SNP IDs algorithm
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
    unsigned i (0);
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
    // reads all the raw-Sample/SNP IDs
    //        vector<string> vecRawSampleIDs;
    //        vector<string> vecRawSNPIDs;
    const unsigned long long FILE_SIZE (0.01*tellFileSize(fileName));
    const unsigned long vecSize(numberSamples*numberSNPs);
    vecGenotypes.reserve(vecSize/4);
    unsigned p1p     (0.015*vecSize);
    //        vecRawSampleIDs.reserve(vecSize);
    //        vecRawSNPIDs.reserve(vecSize);
    unsigned lastp(0);
    unsigned char ctmp(0);
    string stmp;
    unsigned curCharCount(0);

    if(p1p<2)
        p1p=2;
    unsigned lineIndex(0);
    while(getline(infs, dataLine))
    {
        if(++lastp%p1p==0)
            progShow('L', unsigned(infs.tellg()/FILE_SIZE));
        vecStrx=dataLine.split('\t');

        // reading SampleIDs and SNPNames
        if(setSampleIDs.find(vecStrx[1])==setSampleIDs.end())   // new sampleID
        {
            setSampleIDs.insert(vecStrx[1]);
            vecSampleID.push_back(vecStrx[1]);
        }

        if(setSNPIDs.find(vecStrx[0])==setSNPIDs.end())                     // new SNP Name, only 1st sample would go this branch
        {
            setSNPIDs.insert(vecStrx[0]);
            vecSNPName.push_back(vecStrx[0]);
            if(vecStrx[1]!=vecSampleID[0])
            {
                cout << endl;
                E("s", "ERROR: the 1st SNP name of the 2nd sample is not presenting in the 1st sample, bad finalReport file, aborting. X");
                infs.close();
                return -6;
            }
        }
        else if(lineIndex==vecSNPName.size() && lineIndex != numberSNPs)
        {
            infs.close();
            cout << endl;
            I("sis", 3, "SNP IDs read  = ", vecSNPName.size(), "X");
            E("sis", "ERROR: SNP_IDs_read != numberSNPs (", numberSNPs, ") from [Header] section, bad finalReport file, aborting. X");
            return -7;
        }
        else if(vecSNPName[lineIndex % vecSNPName.size()]!=vecStrx[0])      // samples other than 1st
        {
            cout << endl;
            E("s", "ERROR: index out of boundary for nSNPs per animal, please re-check the inputed bad finalReport file, aborting. X");
            E("sssssss", "ERROR: different SNP names/order found, Sample [", vecStrx[1].c_str(),
              "] - SNP [", vecStrx[0].c_str(), "], SNPName should be {", vecSNPName[lineIndex % vecSNPName.size()].c_str(), "}, bad finalReport file, aborting. X");
            infs.close();
            return -6;
        }
        else if(vecSampleID[lineIndex / vecSNPName.size()]!=vecStrx[1])      // samples other than 1st
        {
            cout << endl;
            E("s", "ERROR: index out of boundary for nSNPs per animal, please re-check the inputed bad finalReport file, aborting. X");
            E("sssssss", "ERROR: wrong sample ID found, Sample [", vecStrx[1].c_str(),
              "] - SNP [", vecStrx[0].c_str(), "], sampleID should be {", vecSampleID[lineIndex / vecSNPName.size()].c_str(), "}, bad finalReport file, aborting. X");
            infs.close();
            return -8;
        }

        ++lineIndex;

        // encoding
        stmp=vecStrx[posAllele1AB]+vecStrx[posAllele2AB];
        if(mapAB2uc.find(stmp)==mapAB2uc.end())    // not a valid genotype
        {
            cout << endl;
            E("sssss", "Aborted due to invalid SNP genotype:", stmp.c_str(), "in line", vecStrx[0].c_str(), vecStrx[1].c_str());
            infs.close();
            exit(-1018);
        }
        /* new algorithm for storing genotype vector in the memory
         *
         *  AA    BB    AB/BA NoCall/--
         *  1     2     3     0
         *
         *  each unsigned char holds up to 4 SNPs
         *
         *  unsigned char (5) = 0000 0101 = -- -- AA AA
         *  unsigned char (3) = 0000 0011 = -- -- -- AB
         *  unsigned char (?) = 1110 0100 = AB BB AA --
         */
        ctmp = (ctmp << 2) + mapAB2uc[stmp];   // left shift 2 bits and store current genotype at the lower 2 bits
        //                ctmp+=mapAB2uc[aLine];
        if(++curCharCount==4) // full, plug in and refresh
        {
            vecGenotypes.push_back(ctmp);
            ctmp=curCharCount=0;
        }
    }

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
    if(vecSNPName.size()!=numberSNPs)
    {
        I("sis", 3, "number of SNPNames read = ", vecSNPName.size(), "X");
        E("sis", "ERROR: SNP_Names_read != numberSNPs", numberSNPs, "from [Header] section, bad finalReport file, aborting. X");
    }

    if(setSampleIDs.size()*setSNPIDs.size()!=vecSize)    // check if nSNPs and nSamples consistent with file data body
    {
        E("s", "ERROR: potentially index out of boundary, please re-check if numberSamples and numberSNPs in the header part is consistent with the data body. Bad finalReport file, aborting. X");
        I("si", 1, "[Header] nSNP * nSample =", vecSize);
        I("si", 1, "but tot data entry read =", setSampleIDs.size()*setSNPIDs.size());
        exit (-9);
    }
    set<string>().swap(setSampleIDs);
    set<string>().swap(setSNPIDs);

    for(i=0; i<vecSampleID.size(); ++i)
    {
        for(unsigned j=0; j<vecSampleID[i].size(); ++j)
            if(vecSampleID[i][j]=='#')
                vecSampleID[i][j]='_';
    }
    progClear();
    return 1;
}


//**************************************************************************************************************************

#endif
