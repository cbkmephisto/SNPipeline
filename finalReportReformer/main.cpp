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
//  main.cpp
//  finalReportReformer
//
//  Created by Hailin SU on 10/30/13.
//

#include "../appInfo.h"
#include "../splitablestring.h"
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
    class finalReportReformer : public appInfo
    {
    public:
        finalReportReformer(const int argc, const char * argv[]);
    private:
        unsigned long               startTime;
        string                      inputFileName;
        string                      outputFileName;
        string                      SNPrefFileName;
        string                      headerFileName;
        stringstream                ss;
        bool                        vb;
        // just init once
        vector<string>              vecRefSNPNames;
        // need clear every animal
        map<string, string>         mapSNP2Genotype;
        // methods
        void  refresh();
        bool  parseParameters(const int argc, const char * argv[]);
        void  theMain(const int argc, const char * argv[]);
        void  usage(const string& argv0);
    };

    void finalReportReformer::usage(const string& argv0)
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
        << "        " << argv0 << " " << "-h header -r SNP_Ref -o outputFile -i inputFile" << endl <<endl
        << "options:" << endl
        << "        -r SNP_Ref file name: the expected output SNP_Name template, one SNP name per line, output with the same order" << endl
        << "        -o output  file name" << endl
        << "        -i input   file name: full ill-formatted finalReport file" << endl
        << "        -h header  file name: the expected/corrected header part, including [Header] section, [Data] line and anotation line" << endl
        << endl
        << endl;
    }

    finalReportReformer::finalReportReformer(const int argc, const char * argv[])
    {
        progName    = "finalReportReformer";
        version     = "2014-01-24 15:11 CDT";
        created     = "2013-10-30";
        updated     = "N/A";
        descrip     = "Re-format ill-formatted finalReport files.";
        vb=true;
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
         2016-04-13 18:50CDT    output same SNP order with SNP_Ref file
         */
    }

    bool finalReportReformer::parseParameters(const int argc, const char* argv[])
    {
        //        if(argc < 7)
        //            return false;
        headerFileName = inputFileName = outputFileName = SNPrefFileName = " ";
        string stmp;
        for(unsigned i=1; i<argc; ++i)
        {
            stmp=argv[i];
            if(stmp=="-o") // set the outputFileName
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-o" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                outputFileName=argv[++i];
            }
            else if(stmp=="-r") // set the SNP ref
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-x" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                SNPrefFileName=argv[++i];
            }
            else if(stmp=="-i") // inserts the inputFiles
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                inputFileName=argv[++i];
            }
            else if(stmp=="-h") // the herder file
            {
                if(i+1==argc || argv[i+1][0]=='-')       // if the "-i" goes as the last argv, it is not a full param list.
                {
                    usage(argv[0]);
                    return false;
                }
                headerFileName=argv[++i];
            }
        }
        if(inputFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No input file (-i) specified.");
            return false;
        }
        if(outputFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No output file (-o) specified.");
            return false;
        }
        if(SNPrefFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No SNP ref file (-r) specified.");
            return false;
        }
        if(headerFileName==" ")
        {
            usage(argv[0]);
            E("s", "ERROR: No header file (-h) specified.");
            return false;
        }
        return true;
    }

    void finalReportReformer::theMain(const int argc, const char * argv[])
    {
        startTime = clock();

        // **************** parsing parameters
        if(!parseParameters(argc, argv))
            exit(-777);
        L("s", 1, "Parameters parsed.");
        L("ss", 2, "[Output file Name] : ", outputFileName.c_str());
        L("ss", 2, "[SNPref file Name] : ", SNPrefFileName.c_str());
        L("ss", 2, "[Input  file Name] : ", inputFileName.c_str());
//        set<string>::iterator sit;
//        set<string>::iterator ssit;

        // **************** reading the SNPxref File
        L("s", 1, "Reading the SNP ref File");
        ifstream infs(SNPrefFileName.c_str());
        if(!infs)
        {
            E("sss", "The given SNP_ref file named {", SNPrefFileName.c_str(), "} is invalid. Aborting. X");
            exit(-5);
        }
        splitablestring aLine;
        string stmp;
        while(getline(infs, aLine))  // foreach line of xref
            vecRefSNPNames.push_back(aLine);
        infs.close();   // finish reading the SNP_ref

        // **************** reading the header File
        L("s", 1, "Reading the header file");
        infs.open(headerFileName.c_str());
        if(!infs)
        {
            E("sss", "The given header file named {", headerFileName.c_str(), "} is invalid. Aborting. X");
            exit(-5);
        }
        ofstream ofs (outputFileName.c_str());
        if(!ofs)
        {
            infs.close();
            E("sss", "Couldn't open {", outputFileName.c_str(), "} to write. Aborting. X");
            exit(-5);
        }

        while(getline(infs, aLine))  // foreach line of header, write to output file
            ofs << aLine << endl;
        infs.close();   // finish reading the SNP_ref


        // **************** deal with the input File
        L("s", 1, "Working on input File");
        infs.open(inputFileName.c_str());
        if(!infs)
        {
            E("sss", "The given input file named {", inputFileName.c_str(), "} is invalid. Aborting. X");
            exit(-5);
        }
        // skip to [Data]
        unsigned lineNumber(0);
        while(getline(infs, aLine))  // skip line of header
        {
            ++lineNumber;
            if(aLine=="[Data]")
                break;
        }
        getline(infs, aLine);// skip annotation

        progShow('W', 0);

        const unsigned long long FILE_SIZE(0.01*tellFileSize(inputFileName));
        map<string, string>::iterator mit;
        string curAnm;
        vector<string> sx;
        // 1st line of data
        refresh();  // refresh for the 1st anm
        if(getline(infs, aLine))
        {
            ++lineNumber;
            sx=aLine.split('\t');
            if(sx.size()<4)
            {
                infs.close();
                ofs.close();
                progClear();
                E("si", "Not a 4-col through [Data] section file, line number:", lineNumber);
                exit(-5);
            }

            curAnm=sx[1];
            mapSNP2Genotype[sx[0]]=sx[2]+"\t"+sx[3];    // genotype
//            I("sssss", 1, "created map", sx[0].c_str(), "-> [", mapSNP2Genotype[sx[0]].c_str(), "]");
        }
        // rest of data
        unsigned i;
        while(getline(infs, aLine))
        {
            ++lineNumber;
            sx=aLine.split('\t');
            if(sx.size()<4)
            {
                infs.close();
                ofs.close();
                progClear();
                E("si", "Not a 4-col through [Data] section file, line number:", lineNumber);
                exit(-5);
            }
            if(curAnm!=sx[1])// new animal block
            {
                // write prev anm
//              for(mit=mapSNP2Genotype.begin(); mit!=mapSNP2Genotype.end(); ++mit)
                for(i=0; i<vecRefSNPNames.size(); ++i)
                {
                    stmp=vecRefSNPNames[i];
//                    I("sssss", 1, "writing map", mit->first.c_str(), "-> [", mit->second.c_str(), "]");
                    ofs << stmp << "\t" << curAnm << "\t" << mapSNP2Genotype[stmp] << endl;
                }
                // go on
                curAnm=sx[1];
                progShow('W', unsigned(infs.tellg()/FILE_SIZE));
                refresh();  // refresh for the next anm
            }
            mapSNP2Genotype[sx[0]]=sx[2]+"\t"+sx[3];    // genotype
        }

        // write out the last anm
//        for(mit=mapSNP2Genotype.begin(); mit!=mapSNP2Genotype.end(); ++mit)
//            ofs << mit->first << "\t" << curAnm << "\t" << mit->second << endl;

        for(i=0; i<vecRefSNPNames.size(); ++i)
        {
            stmp=vecRefSNPNames[i];
            ofs << stmp << "\t" << curAnm << "\t" << mapSNP2Genotype[stmp] << endl;
        }
        map<string, string>().swap(mapSNP2Genotype);
        progClear();
        infs.close();   // finish reading the SNPxref
        L("sis", 1, "Done", (clock()-startTime)/CLOCKS_PER_SEC, "seconds." );
        // deal with
    }

    void finalReportReformer::refresh()
    {
        map<string, string>().swap(mapSNP2Genotype);
        for(unsigned i=0; i<vecRefSNPNames.size(); ++i)
            mapSNP2Genotype[vecRefSNPNames[i]]="-\t-";
    }
}

int main(int argc, const char * argv[])
{
    libcbk::finalReportReformer(argc, argv);
    return 0;
}
