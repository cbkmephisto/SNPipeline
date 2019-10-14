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
//  SNPipeline.mergingAB.cpp
//  libHailins
//
//  Created by Hailin SU on 11/22/2013.
//
//  SNPipeline - step 2 - merging AB genotype files

#include "appInfo.h"
#include "splitablestring.h"
#include <vector>
#include <map>
#include <set>
#include <iomanip>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)
#define L(...)              if(vb) lg.generalLog("I", __VA_ARGS__)

// comment out the next line to undefine D_EBUG and disable all the outputs with L()
/*#define D_EBUG

#ifdef D_EBUG
    #define L(...)              lg.generalLog("I", __VA_ARGS__)
#else
    #define L(...)
#endif
*/
using namespace std;
namespace libcbk
{
    class SNPipeline_mergingAB : public appInfo
    {
    public:
        SNPipeline_mergingAB(const int argc, const char * argv[]);
    private:
        string                      outputFileName;
        string                      SNPxrefFileName;
        string                      genotypeTmp;
        stringstream                ss;
        bool                        vb;
        float                       crLimit;
        // just init once
        set<string>                 setInputFileNames;
        map<string, string>         mapSNPxref;
        map<string, unsigned long>  mapSNPName2Pos;
        set<string>                 setSampleIDs;
        set<string>                 setDuplicatedSampleIDs;
        vector<string>              vecDuplicatedSampleIDs;
        vector< vector<string> >    vecDuplicatedABGenotypes;
        vector<string>              sx;
        vector<string>              vecDiscLQ_id;
        vector<string>              vecDiscLQ_from;
        vector<float>               vecDiscLQ_cr;
        // need clear for each input file
        vector<string>              vecOptABGenotype;
        vector<string>              vecSNP_IDsLine;

        // methods
        bool  parseParameters(const int argc, const char * argv[]);
        void  theMain(const int argc, const char * argv[]);
        void  usage(const string& argv0);
        float callingRate(const vector<string>& vec, const unsigned len);
        float similarity(const vector<string>& vec1, const vector<string>& vec2);
    };

    void SNPipeline_mergingAB::usage(const string& argv0)
    {
        cout << endl
        << "     version: " << version << endl
        << "####################################################################" << endl
        << "#######" << endl
        << "###                 " << progName << endl
        << "####" << endl
        << "############" << endl
        << "##                              by hailins@iastate.edu" << endl
        << "#####                           created " << created << endl
        << "####                            updated " << updated << endl
        << "#######" << endl
        << "####################################################################" << endl
        << endl << "    " << descrip << endl
        << endl
        << "usage:" << endl
        << "        1) " << argv0 << " " << "-x SNPxref -o outputFile -l list_file [-v] [-r 0.95]" << endl
        << "        2) " << argv0 << " " << "-x SNPxref -o outputFile -i AB-genotypeFile1 [AB-genotypeFile2 *] [-v] [-r 0.95]" << endl
        << "        3) " << argv0 << " " << "-x SNPxref -o outputFile -i AB-genotypeFile1 [AB-genotypeFile2 *] -l list_file [-v] [-r 0.95]" << endl <<endl
        << "options:" << endl
        << "        -x SNPxref file name"<< endl
        << "        -o output file name" << endl
        << "        -i input file name (list)" << endl
        << "        -l input list file containing input files, seperated by new line" << endl
        << "        -v verbose output" << endl
        << "        -r 0~1, calling rate limitation: only write out samples with calling rate greater than/o/equal to the given value." << endl
        << "           Not specifying this value is the same as specifying 0." << endl
        << endl
        << "The order of [-r], [-v], -x block, -i block, -l block and -o block doesn't matter." << endl
        << endl;
    }

    SNPipeline_mergingAB::SNPipeline_mergingAB(const int argc, const char * argv[])
    {
        progName    = "SNPipeline.mergingAB";
        version     = "2014-08-07 15:40 CDT";
        created     = "2013-11-22";
        updated     = "2014-08-07";
        descrip     = "SNPipeline - setp 2. Creates merged-AB-genotype-file from AB-genotype file(s).";
        //        displayInfo();
        //        lg.initLogger();
        theMain(argc, argv);
        //        lg.finalizeLogger();
        /*
         updates
         2013-11-27 15:11 CDT   local vector<string> vecPosX get rid of slow calling vec[map[vec[j]]]
                                    duplicated samples were compared using calling rate and templates
         2014-01-24 15:09 CDT   - fixed potential bugs in recognizing parameters and SampleIDs
                                - SNPxref can be space/tab delimited
         2014-02-06 11:00 CDT   calling rate limitation and merging duplicated samples with genotypic similarity info
         2014-02-07 12:07 CDT   filtering calling rate and genotypic similarity, write the bad outcomes into files
         2014-02-21 11:25 CDT   output duplicated samples with identical call rate to .hiSim file
         2014-08-07 15:39 CDT   output a new file for ID - path pairs.
         */
    }

    bool SNPipeline_mergingAB::parseParameters(const int argc, const char* argv[])
    {
//        if(argc < 7)
//            return false;
        outputFileName = SNPxrefFileName = " ";
        crLimit=0.0;
        string stmp;
        for(unsigned i=1; i<argc; ++i)
        {
            stmp=argv[i];
            if(stmp=="-l")      // read the list
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-l" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                stmp=argv[++i];     // stores the list-fileName and open it for reading
                ifstream infs(stmp.c_str());
                if(!infs)
                {
                    E("sss", "The given listFileName {", argv[i], "} is invalid. Aborting. X");
                    exit(-1);
                }
                while(getline(infs, stmp))
                    setInputFileNames.insert(stmp);
                infs.close();
            }
            else if(stmp=="-v") // set verbose output
                vb=true;
            else if(stmp=="-o") // set the outputFileName
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-o" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                outputFileName=argv[++i];
            }
            else if(stmp=="-x") // set the SNPxref
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-x" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                SNPxrefFileName=argv[++i];
            }
            else if(stmp=="-r") // calling rate limitation
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-r" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                string cr_str (argv[++i]);
                stringstream ss;
                ss << cr_str;
                ss >> crLimit;
                ss.clear();
            }
            else if(stmp=="-i") // inserts the inputFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                for(++i; i<argc && argv[i][0]!='-'; ++i)
                {
                    setInputFileNames.insert(argv[i]);
                }
                if(i!=argc)         // if not the end of param list, put the index back 1 position.
                    --i;
            }
        }
        if(setInputFileNames.size()<1)
        {
            usage(argv[0]);
            E("s", "ERROR: No input file(s) (-i/-l) specified.");
            return false;
        }
        if(outputFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No output file (-o) specified.");
            return false;
        }
        if(SNPxrefFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No SNPxref file (-x) specified.");
            return false;
        }
        return true;
    }

    float SNPipeline_mergingAB::callingRate(const vector<string> &vec, const unsigned len)
    {
        static unsigned i;
        static float ret;
        ret=vec.size();
        for(i=0; i<vec.size(); ++i)
            if(vec[i]=="--")
                --ret;
        return ret/(float)(len);
    }

    float SNPipeline_mergingAB::similarity(const vector<string>& vec1, const vector<string>& vec2)
    {
        static unsigned i, len, clen;
        static float ret;
        ret=0;
        clen=len=(unsigned)vec1.size();
        for(i=0; i<clen; ++i)
            if (vec1[i]=="--" || vec2[i]=="--")    // noCall
                --len;
            else if(vec1[i]==vec2[i])       // called, and eq
                ++ret;
        return ret/(float)(len);
    }

    void SNPipeline_mergingAB::theMain(const int argc, const char * argv[])
    {
        // **************** parsing parameters
        if(!parseParameters(argc, argv))
            exit(-777);
        L("s", 1, "Parameters parsed.");
        L("ss", 2, "[Output file Name ] : ", outputFileName.c_str());
        L("ss", 2, "[SNPxref file Name] : ", SNPxrefFileName.c_str());
        L("s",  2, "[Input files list ] :");
        set<string>::iterator sit;
        set<string>::iterator ssit;
        for(sit=setInputFileNames.begin(); sit!=setInputFileNames.end(); ++sit)
            L("s", 3, (*sit).c_str());

        // **************** reading the SNPxref File
        L("s", 1, "Reading the SNPxref File");
        ifstream infs(SNPxrefFileName.c_str());
        if(!infs)
        {
            E("sss", "The given SNPxref file name {", SNPxrefFileName.c_str(), "} is invalid. Aborting. X");
            exit(-5);
        }
        splitablestring aLine;
        string stmp;
        unsigned currentPos(0);
        while(getline(infs, aLine))  // foreach line of xref
        {
//            aLine.erase(remove(aLine.begin(), aLine.end(), '\r'), aLine.end()); // get rid of the ^M char
            ss << aLine;            // counts the columns
            while(ss>>stmp)
                sx.push_back(stmp);
            ss.clear();
            if(sx.size()!=2)
            {
                // wrong xref file, skipping
                E("s", "ERROR: illegal SNPxref file read (not a 2-column-whole-through space/tab-delimited file). Aborting. X");
                cerr << "              line content: " << aLine << endl;
                map<string, string>().swap(mapSNPxref);
                infs.close();
                exit(-123);
            }
            mapSNPxref[sx[0]]=sx[1];
            if(mapSNPName2Pos.find(sx[1])==mapSNPName2Pos.end())    // push new SNP name
                mapSNPName2Pos[sx[1]]=currentPos++;
            vector<string>().swap(sx);
        }
        infs.close();   // finish reading the SNPxref
        vecOptABGenotype.resize(mapSNPName2Pos.size());   // determine the size of genotype line

        // **************** read the header of all the input files to check if SNPxref contains all the SNP names in all the input AB-genotype file and fills all the Sample IDs to check duplications
        L("s", 1, "Reading input files to check SNP_IDs and duplicated sample IDs");
        unsigned i(0), j(0);
        bool foundHeader(false);
        for(sit=setInputFileNames.begin(); sit!=setInputFileNames.end(); ++sit)
        {
            L("ss", 2, "reading", (*sit).c_str());
            infs.open((*sit).c_str());
            if(!infs)
            {
                E("sss", "The given input file name {", (*sit).c_str(), "} is invalid. Aborting. X");
                exit(-1*i-100);
            }
            unsigned long FILE_SIZE(0.01*tellFileSize(*sit));
            if(FILE_SIZE==0) FILE_SIZE=100;
            foundHeader=false;
            if(getline(infs, aLine) && aLine[0]=='#')
            {
                sx=aLine.split('#');
                if(sx[0]=="SNP_IDs")
                    foundHeader=true;
            }
            if(!foundHeader)    // 10 lines not founding header info
            {
                infs.close();
                E("sss", "The given input AB-genotype file {", (*sit).c_str(), "} is invalid. Aborting. X");
                exit(-1*i-200);
            }
            vector<string> vecInvalidSNPname;
            for(j=1; j<sx.size(); ++j)  // check the SNP_Names. sx[0]="", sx[1]="SNP_IDs"
            {
                if(mapSNPxref.find(sx[j])==mapSNPxref.end())    // the SNP_Name is not in the SNPxref list
                    vecInvalidSNPname.push_back(sx[j]);
            }
            if(vecInvalidSNPname.size()!=0)
            {
                infs.close();
                E("sss", "The given input AB-genotype file {", (*sit).c_str(), "} contains illegal SNP Name(s):");
                for(i=0; i<vecInvalidSNPname.size(); ++i)
                    cout << "\t\t" << vecInvalidSNPname[i] << endl;
                E("s", "Aborting. X");
                exit(-1*i-1000000);
            }
            // fills vec/set-SampleIDs
            unsigned linesRead(0);
            ss.str("");
            ss.clear();
            while(getline(infs, aLine))
            {
                if(aLine[0]=='#')   // path info line
                    continue;
                ss << aLine;
                ss >> stmp;    // get the 1st column: the Sample ID
                ss.str("");
                ss.clear();
                if(setSampleIDs.find(stmp)==setSampleIDs.end())    // new/non-duplicated
                    setSampleIDs.insert(stmp);
                else                                                // duplicated found
                    setDuplicatedSampleIDs.insert(stmp);
//                getline(infs,aLine);    // rest of the line
                if(++linesRead%100==1)
                    progShow('R', unsigned(infs.tellg()/FILE_SIZE));
            }
            infs.close();
            progClear();
        }
        if(setDuplicatedSampleIDs.size()!=0)
        {
            E("sis", "WARNING:", setDuplicatedSampleIDs.size(), "duplicated sample ID(s) found, merge later.");
            // for each duplicated Sample ID,
//            for(sit=setDuplicatedSampleIDs.begin(); sit!=setDuplicatedSampleIDs.end(); ++sit)
//                cout << (*sit) << " ";
//            cout << endl;
        }

        // **************** start merging
        L("s", 1, "Start merging");
        ofstream ofs (outputFileName.c_str());
        if(!ofs)
        {
            E("sss", "Initializing output file {", outputFileName.c_str(), "} failed. Aborting. X");
            exit(-10);
        }
        string fnID2Path (outputFileName + ".ID2Path");
        ofstream ofsID2Path (fnID2Path.c_str());
        ofsID2Path << "#STATUS:" << endl;
        ofsID2Path << "# -good- : this entry was in the result file, no duplicated & high call rate" << endl;
        ofsID2Path << "# low_CR : this entry was discarded due to low call rate" << endl;
        ofsID2Path << "# dupQFC : this entry had qualified call rate but was duplicated with other entry(ies)" << endl;
        ofsID2Path << "# dupLCR : this entry was discarded due to low call rate, but was duplicated with other entry(ies)" << endl<<endl;
        ofsID2Path << "#ID\tstatus\tpath" << endl;
        
        L("s", 2, "Writing out the headerLine (#SNP_IDs#)");
        ofs << "#SNP_IDs";
        vector<string>().swap(sx);
        sx.resize(mapSNPName2Pos.size());
        map<string, unsigned long>::iterator itMapSNPName2Pos;
        for(itMapSNPName2Pos=mapSNPName2Pos.begin(); itMapSNPName2Pos!=mapSNPName2Pos.end(); ++itMapSNPName2Pos)
            sx[itMapSNPName2Pos->second]=itMapSNPName2Pos->first;
        for(i=0; i<sx.size(); ++i)
            ofs << "#" << sx[i];
        ofs << endl;
        vector<string>().swap(sx);
        float cr;

        L("s", 2, "Entering inputFileName loop");
        string curPath;
        vector<string> vecDupPath;
        map<string, unsigned> mapPath2numSNPs;
        unsigned numberSNPs;
        for(sit=setInputFileNames.begin(); sit!=setInputFileNames.end(); ++sit)
        {
            L("ss", 3, "Processing", (*sit).c_str());
            infs.open((*sit).c_str());
            unsigned long FILE_SIZE(0.01*tellFileSize(*sit));
            if(FILE_SIZE==0) FILE_SIZE=100;
            // clear
            vector<string>().swap(vecSNP_IDsLine);
            // 1st line: SNP_IDs
            getline(infs, aLine);
            vecSNP_IDsLine=aLine.split('#');    //  vecSNP_IDsLine[0]=="SNP_IDs", [1]==0th SNP name, etc
            vector<long> vecPosX;
            vecPosX.resize(vecSNP_IDsLine.size()-1);
            numberSNPs=unsigned(vecPosX.size());    // numberSNPs is input-AB-genotype file specific
            for(j=1; j<vecSNP_IDsLine.size(); ++j) // rename the SNP_IDs and store the output pos for genotypes
            {
                vecSNP_IDsLine[j]=mapSNPxref[vecSNP_IDsLine[j]];
                vecPosX[j-1]=mapSNPName2Pos[vecSNP_IDsLine[j]];
            }
            // throughout the file
            unsigned linesRead(0);
            while(getline(infs, aLine))
            {
                if(aLine[0]=='#')   // the path line
                {
                    // sx=aLine.split('#');
                    curPath = aLine;
                    continue;
                }
                sx=aLine.split(' ');    //          sx[0]==Sample_ID, [1]==0th SNP AB genotype, etc
                for(j=0; j<vecOptABGenotype.size(); ++j)    // empty & fill the genotype
                    vecOptABGenotype[j]="--";
                // populates the genotype
                for(j=1; j<sx.size(); ++j)
                    vecOptABGenotype[vecPosX[j-1]] = sx[j];
                // write out
                if(setDuplicatedSampleIDs.find(sx[0])==setDuplicatedSampleIDs.end()) // not a duplicated one
                {
                    // test calling rate limitation
                    if(crLimit==0.0 || (cr=callingRate(vecOptABGenotype, numberSNPs)) >= crLimit)    // qualified output
                    {
                        // construct the genotype line string
                        genotypeTmp="";
                        for(j=0; j<vecOptABGenotype.size(); ++j)
                        {
                            genotypeTmp+=" ";
                            genotypeTmp+=vecOptABGenotype[j];
                        }
                        ofs << setfill(' ') << left << setw(32) << sx[0] << genotypeTmp << endl;
                        ofsID2Path << sx[0] << "\t-good-\t" << curPath << endl;
                    }
                    else   // discard
                    {
                        vecDiscLQ_id.push_back(sx[0]);
                        vecDiscLQ_from.push_back(curPath);
                        vecDiscLQ_cr.push_back(cr);
                        ofsID2Path << sx[0] << "\tlow_CR\t" << curPath << endl;
//                        E("ssf", "        discarding low-calling-rate sample", sx[0].c_str(), cr);
                    }
                }
                else        // duplicated, discard low CR ones, keep high CR to merge later
                {
                    if(crLimit==0.0 || (cr=callingRate(vecOptABGenotype, numberSNPs)) >= crLimit)    // qualified duplication
                    {
                        vecDuplicatedSampleIDs.push_back(sx[0]);
                        vecDuplicatedABGenotypes.push_back(vecOptABGenotype);//(vector<string>(vecOptABGenotype));
                        vecDupPath.push_back(curPath);
                        mapPath2numSNPs[curPath]=numberSNPs;
                        ofsID2Path << sx[0] << "\tdupQFC\t" << curPath << endl;
                    }
                    else   // discard
                    {
                        vecDiscLQ_id.push_back(sx[0]);
                        vecDiscLQ_from.push_back("DUP_LQ:" + curPath);
                        vecDiscLQ_cr.push_back(cr);
                        ofsID2Path << sx[0] << "\tdupLCR\t" << curPath << endl;
                        //                        E("ssf", "        discarding low-calling-rate sample", sx[0].c_str(), cr);
                    }
                }
                if(++linesRead%100==1)
                {
//                    vector<string>().swap(sx);
                    progShow('P', unsigned(infs.tellg()/FILE_SIZE));
                }
            }   // end of getline loop
            infs.close();
            progClear();
        }   // end of foreach input fileName
        ofsID2Path.close();

        // deals with duplicated samples
        // 0 these are all qualified calling rate duplicated animals
        // 1 find the genotype with highest call rate as template
        // 2 fills missing genotypes of template if possible
        // 3 if anyone of current duplicated animalID got a low sim, discard to lowSim

        const float simThreshold (0.9);
        unsigned vi, vj;    // for vecDuplicatedIndex[i] and j
        if(setDuplicatedSampleIDs.size()!=0)
        {
            ofstream ofsLowSim, ofsHiSim;
            bool blowSim(false), bhiSim(false);
            vector<unsigned> vecDuplicatedIndex;
            float sim, simMin(2.0), simMax(0.0);
            L("s", 2, "Merging duplicated samples");
            // for each duplicated Sample ID,
            for(ssit=setDuplicatedSampleIDs.begin(); ssit!=setDuplicatedSampleIDs.end(); ++ssit)
            {
                // #### clear the vecIndex
                vector<unsigned>().swap(vecDuplicatedIndex);

                L("s", 3, ssit->c_str());
                vecOptABGenotype[0]="ST";   // mark as new loop
                // #### 1st loop, determine the template, and store the genotypes for similarity check
                for(i=0; i<vecDuplicatedSampleIDs.size(); ++i)
                {
                    if(vecDuplicatedSampleIDs[i]!=(*ssit))   // skipping others, just deal with the current iteration
                        continue;
                    cr=callingRate(vecDuplicatedABGenotypes[i], mapPath2numSNPs[vecDupPath[i]]);
//                    I("sfss", 4, "[", cr, "]\tfrom #", vecPath[i].c_str());
                    if(vecOptABGenotype[0]=="ST")   // 1st, just assign
                        vecOptABGenotype=vecDuplicatedABGenotypes[i];
                    else                            // others, compare the calling rate, get the larger as template
                    {
                        if(callingRate(vecOptABGenotype, mapPath2numSNPs[vecDupPath[i]])<cr)
                            vecOptABGenotype=vecDuplicatedABGenotypes[i];
                    }
                    vecDuplicatedIndex.push_back(i); // store the genotype index for similarity check
                }

                // #### 2nd loop, filling out the gaps in the template if possible
                for(i=0; i<vecDuplicatedSampleIDs.size(); ++i)
                {
                    if(vecDuplicatedSampleIDs[i]!=(*ssit))   // skipping others, just deal with the current iteration
                        continue;
                    for(j=0; j<vecOptABGenotype.size(); ++j)    // fill from the stored genotype
                    {
                        if(vecOptABGenotype[j]=="--" && vecDuplicatedABGenotypes[i][j]!="--")
                            vecOptABGenotype[j] = vecDuplicatedABGenotypes[i][j];
                    }
                }

                // #### quality control filter : similarity test.\
                            Failing in this test suggests that the duplicated \
                                animals may not be pointing to the same inidividual.
                // for each genotype vector of current duplicated sample, get the min/manx similarity
                simMin=2;
                simMax=0;
                for(i=0; i<vecDuplicatedIndex.size(); ++i)
                    for(j=0; j<=i; ++j)
                    {
                        vi=vecDuplicatedIndex[i];
                        vj=vecDuplicatedIndex[j];
                        aLine=vecDuplicatedSampleIDs[vi];
                        stmp=vecDuplicatedSampleIDs[vj];
                        if( aLine==(*ssit) && aLine==stmp )// skip if not current animal
                        {
                           if( (sim=similarity(vecDuplicatedABGenotypes[vi], vecDuplicatedABGenotypes[vj]) )
                              <simMin)
                               simMin=sim;
                           else if(sim>simMax)
                               simMax=sim;
                        }
                    }

                if(simMin<simThreshold)  // low SIM, write out
                {
                    if(!ofsLowSim)
                    {
                        ofsLowSim.open((outputFileName+".lowSim").c_str());
                        blowSim=true;
                    }
                    static bool headed;
                    headed=false;
                    for(i=0; i<vecDuplicatedIndex.size(); ++i)
                    {
                        vi=vecDuplicatedIndex[i];
                        aLine=vecDuplicatedSampleIDs[vi];
                        if (aLine!=(*ssit))  // skip if not current animal
                            continue;
                        if(!headed)
                        {
                            ofsLowSim << aLine << endl;
                            headed=true;
                        }
                        for(j=0; j<i; ++j)
                        {
                            vj=vecDuplicatedIndex[j];
                            stmp=vecDuplicatedSampleIDs[vj];
                            if (stmp!=aLine)  // skip if not current animal
                                continue;
                            ofsLowSim << "\t" << similarity(vecDuplicatedABGenotypes[vi], vecDuplicatedABGenotypes[vj]);
                        }
                        curPath=vecDupPath[vi];
                        ofsLowSim << "\t" << 1 << "\tCR "
                        << callingRate(vecDuplicatedABGenotypes[vi], mapPath2numSNPs[curPath])
                        << " from # " << curPath << endl;
                    }
                }
                else
//                if(simMin>=simThreshold && (crLimit==0.0 || cr >= crLimit) )    // qualified output
                {
                    // construct the genotype line string
                    genotypeTmp="";
                    for(j=0; j<vecOptABGenotype.size(); ++j)
                    {
                        genotypeTmp+=" ";
                        genotypeTmp+=vecOptABGenotype[j];
                    }
                    ofs << setfill(' ') << left << setw(32) << (*ssit) << genotypeTmp << endl;
                }
                /*                else if(simMin>=simThreshold)  // discard low calling rate & sim
                 {
                 vecDiscLQ_id.push_back(*ssit);
                 vecDiscLQ_from.push_back("multiple_entry");
                 vecDiscLQ_cr.push_back(cr);
                 //                    E("ssf", "        discarding low-calling-rate duplicated sample after merging:", (*ssit).c_str(), cr);
                 }
                 */

                if(simMax==1)  // high SIM, write out
                {
                    if(!ofsHiSim)
                    {
                        ofsHiSim.open((outputFileName+".idSim").c_str());
                        bhiSim=true;
                    }
                    static bool headed;
                    headed=false;
                    for(i=0; i<vecDuplicatedIndex.size(); ++i)
                    {
                        vi=vecDuplicatedIndex[i];
                        aLine=vecDuplicatedSampleIDs[vi];
                        if (aLine!=(*ssit))  // skip if not current animal
                            continue;
                        if(!headed)
                        {
                            ofsHiSim << aLine << endl;
                            headed=true;
                        }
                        for(j=0; j<i; ++j)
                        {
                            vj=vecDuplicatedIndex[j];
                            stmp=vecDuplicatedSampleIDs[vj];
                            if (stmp!=aLine)  // skip if not current animal
                                continue;
                            ofsHiSim << "\t" << similarity(vecDuplicatedABGenotypes[vi], vecDuplicatedABGenotypes[vj]);
                        }
                        curPath=vecDupPath[vi];
                        ofsHiSim << "\t" << 1 << "\tCR "
                        << callingRate(vecDuplicatedABGenotypes[vi], mapPath2numSNPs[curPath])
                        << " from # " << curPath << endl;
                    }
                }
            }       // end for each setDuplicatedSampleIDs
            if(blowSim)
            {
                E("sss", "Animals with LOW GENOTYPIC SIMILARITY were writen into [", (outputFileName+".lowSim").c_str(), "]");
                ofsLowSim.close();
            }
            if(bhiSim)
            {
                E("sss", "Animals with IDENTICAL GENOTYPEs were writen into [", (outputFileName+".idSim").c_str(), "]");
                ofsHiSim.close();
            }
        }           // end if(setDuplicatedSampleIDs.size()!=0)
        
        // deals with low calling rate samples
        if(vecDiscLQ_id.size()!=0)
        {
            ofstream ofsLQ((outputFileName+".lowCR").c_str());
            ofsLQ << left << setfill(' ')
            << setw(24) << "# AnimalID"
            << setw(16) << "Calling_Rate" << "From" << endl;
            for(i=0; i<vecDiscLQ_id.size(); ++i)
                ofsLQ << left << setfill(' ')
                << setw(24) << vecDiscLQ_id[i]
                << setw(16) << vecDiscLQ_cr[i]
                << vecDiscLQ_from[i] << endl;
            ofsLQ.close();
            E("sss", "Animals with LOW CALLING RATE were writen into [", (outputFileName+".lowCR").c_str(), "]");
        }

        ofs.close();
        L("s", 1, "Done.");
    }
}

int main(int argc, const char * argv[])
{
    libcbk::SNPipeline_mergingAB(argc, argv);
    return 0;
}
