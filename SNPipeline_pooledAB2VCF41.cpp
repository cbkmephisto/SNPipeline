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

//  SNPipeline_pooledAB2VCF41.cpp
//  SNPipeline
//
//  Created by Hailin SU on 2/11/14.
//
//  SNPipeline - step 4 - convert merged/unmerged AB genotype files to VCF4.1 format

#include "appInfo.h"
#include "SNP.h"
#include "splitablestring.h"
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <iomanip>

#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)
#define L(...)              if(vb) lg.generalLog("I", __VA_ARGS__)

using namespace std;

namespace libcbk
{

    class SNPipeline_pooledAB2VCF41 : public appInfo
    {
    public:
        SNPipeline_pooledAB2VCF41(const int argc, const char * argv[]);
    private:

        vector<SNP> vecSNPinMap;

        map<string, string> mapSNPName2Chr;
        map<string, string> mapSNPName2Pos;

        vector<string> vecSNPNames;
        vector<string> vecSampleIDs;

        map<string, unsigned>       mapAB2uc;
        vector<string>              vecUC2STR;

        map<string, string>         mapSNP2Line;
        vector<unsigned char> vecGenotypes;
        bool vb;
        // methods
        void  theMain(const int argc, const char * argv[]);
        void  usage(const string& argv0);
    };

    void SNPipeline_pooledAB2VCF41::usage(const string& argv0)
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
        << argv0 << " mapFile AB-GenotypeFile" << endl << endl
        << "mapFile"<<endl
        << "  3-column file, SNPName chr pos, whitespace/tab delimited" <<endl
        << "  NO HEADER LINE, or headerline starts with '#'. ANY '#' line would be ignored while making the map." << endl << endl
        << "AB-GenotypeFile" << endl
        << "  merged/unmerged ab-genotype file, 1st line containing 'SNP_IDs' and SNPNames delimited by #. Any other lines starting with # would be discarded."<< endl << endl
        << "Discarding SNPs if" << endl
        << "  - SNP not in map file" << endl
        << "  - SNP genotype is all missing in the current file" << endl
        << "  - SNPs sharing the same position, only 1 would be kept"
        << endl     << endl;
    }

    SNPipeline_pooledAB2VCF41::SNPipeline_pooledAB2VCF41(const int argc, const char * argv[])
    {
        progName    = "SNPipeline.pooledAB2VCF41";
        version     = "2014-02-07 11:33 CDT";
        created     = "2014-02-22";
        updated     = "N/A";
        descrip     = "SNPipeline - setp 4. Convert merged/unmerged AB genotype files to VCF4.1 format.";
        //        displayInfo();
        //        lg.initLogger();
        theMain(argc, argv);
        //        lg.finalizeLogger();
        /*
         updates
         2014-02-11 11:33 AM init
         */
    }

    void SNPipeline_pooledAB2VCF41::theMain(const int argc, const char * argv[])
    {
        // **************** parsing parameters
        if(argc!=3)
        {
            usage(argv[0]);
            exit (-1);
        }
        string mapFileName (argv[1]);
        string abgFileName (argv[2]);
        vb=true;

        // **************** reading the map File
        L("s", 1, "Reading the map File");
        ifstream infs(mapFileName.c_str());
        if(!infs)
        {
            E("sss", "The given map file name {", mapFileName.c_str(), "} is invalid. Aborting. X");
            exit(-2);
        }
        splitablestring aLine;
        string stmp;
        stringstream ss;
        unsigned long ltmp;
        unsigned utmp;
        vector<string> sx;
        unsigned long FILE_SIZE(0.01*tellFileSize(mapFileName));
        if(FILE_SIZE<100)
            FILE_SIZE=100;
        unsigned linesRead(0);
        while(getline(infs, aLine))  // foreach line of mapFile
        {
//            ss << aLine;            // counts the columns
//            while(ss>>stmp)
//                sx.push_back(stmp);
//            ss.clear();
//            I("s", 1, aLine.c_str());
//            L("s", 4, aLine.c_str());
            sx=aLine.split();
            if(sx.size()!=3)
            {
                // wrong map file
                cout << endl;
                E("s", "ERROR: illegal map file read (not a 3-column-whole-through space/tab-delimited file). Aborting. X");
                cerr << "              line content: " << aLine << endl;
                infs.close();
                exit(-123);
            }
            // skip # line, for header
            if(sx[0][0]=='#')
            {
//                L("s", 5, "skipping");;
                continue;
            }
            ss.str("");
            ss << sx[1];
            ss >> utmp;     // chr
            ss.clear();
            ss << sx[2];
            ss >> ltmp;     // pos
            ss.clear();

            // skip chr0, pos0 pool SNP
            if(utmp==0 || 0==ltmp)
            {
//                L("s", 5, "skipping");;
                continue;
            }

            mapSNPName2Chr[sx[0]]=sx[1];
            mapSNPName2Pos[sx[0]]=sx[2];
            vecSNPinMap.push_back(SNP(sx[0], utmp, ltmp));
//            vector<string>().swap(sx);
            if(++linesRead%1000==1)
                progShow('M', unsigned(infs.tellg()/FILE_SIZE));
        }
        progClear();
        infs.close();   // finish reading the mapFile
        L("is", 2, mapSNPName2Chr.size(), "SNPs read from the map.");
        // **************** read the genotype
        L("s", 1, "Reading genotype file");
        infs.open(abgFileName.c_str());
        if(!infs)
        {
            E("sss", "The given ab-genotype file name {", abgFileName.c_str(), "} is invalid. Aborting. X");
            exit(-100);
        }
        unsigned i(0), j(0);
        FILE_SIZE = 0.01*tellFileSize(abgFileName);
        if(FILE_SIZE==0) FILE_SIZE=100;


        // 1st line, vecSNPNames
        if(getline(infs, aLine))
            vecSNPNames=aLine.split('#');   // vecSNPNames[0]=="SNP_IDs", remember 1st SNPName is vecSNPNames[1]

        L("is", 2, vecSNPNames.size()-1, "SNPs read from the ab-genotype file.");

        progShow('L', 0);

        // check SNPName in map
        set<string> setDiscardedSNP;
        for(i=1; i<vecSNPNames.size(); ++i)
            if(mapSNPName2Chr.find(vecSNPNames[i])==mapSNPName2Chr.end())
                setDiscardedSNP.insert(vecSNPNames[i]);

        // init static map/vector
        if(mapAB2uc.size()==0)
        {
            mapAB2uc["--"]=0;
            mapAB2uc["AA"]=1;
            mapAB2uc["BB"]=2;
            mapAB2uc["AB"]=3;
        }
        if(vecUC2STR.size()==0)
        {
            vecUC2STR.push_back("./.");// 0, --
            vecUC2STR.push_back("0/0");// 1, AA
            vecUC2STR.push_back("1/1");// 2, BB
            vecUC2STR.push_back("0/1");// 3, AB
        }
        // rest of the lines of ab-genotype file, fills the genotype vector and sampleID vector
        linesRead=0;
        while(getline(infs, aLine))
        {
            if(++linesRead%100==1)
            {
                progShow('L', unsigned(infs.tellg()/FILE_SIZE));
                ss.str("");
                ss.clear();
            }
            if(aLine[0]=='#')   // path info line
                continue;
            ss << aLine;
            ss >> stmp;     // get the 1st column: the Sample ID
            vecSampleIDs.push_back(stmp);

            while(ss>>stmp) // other columns: genotypes for current sample
                vecGenotypes.push_back(mapAB2uc[stmp]);
            ss.clear();
        }
        infs.close();
        progClear();
        const unsigned long numberSNPs      (vecSNPNames.size()-1);
        const unsigned long numberSamples   (vecSampleIDs.size());

        // **************** creating SNP_ID -> output line map
        L("s", 2, "creating SNP_ID -> output line map");
        // foreach SNP as an output line
        unsigned tn (numberSNPs*.01);
        if(tn<100)
            tn=100;
        unsigned gtSum;
        set<string> setDiscardedSNPMS;
        unsigned char ctmp(0);
        for(i=0; i<numberSNPs; ++i)
        {
            if(i%tn==0)
                progShow('C', i/tn);

            stmp=vecSNPNames[i+1];  // vecSNPNames[0]=="SNP_IDs"
            // check if SNP NOT in map? continue
            if(setDiscardedSNP.find(stmp)!=setDiscardedSNP.end())
            {
                //                cout << endl << "WARNING: discarded SNP_ID not in map file: " << stmp << endl;
                continue;
            }
            // check if all SNPs are missing? continue: print
            gtSum=0;
            for(j=0; j<numberSamples; ++j)
                //                  base         + offset
                gtSum+=vecGenotypes[j*numberSNPs + i     ]; // 0 for missing
            if(gtSum==0)
            {
                //                cout << endl << "WARNING: discarded SNP genotype all missing in ab-genotype file: " << stmp << endl;
                setDiscardedSNPMS.insert(stmp);
                continue;
            }
            // print meta info #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
//            aLine = mapSNPName2Chr[stmp]        // CHROM
//            + "\t" + mapSNPName2Pos[stmp]       // POS
//            + "\t" + stmp                       // ID
//            + "\tG\tA\t.\tPASS\t.\tGT";         // REF ALT QUAL. PASS INFO. FORMAT
            // foreach animal as genotype output column
            aLine=string("");
            ctmp=0;
            for(j=0; j<numberSamples; ++j)
            {
                ctmp = (ctmp << 2) + vecGenotypes[j*numberSNPs + i];   // left shift 2 bits and store current genotype at the lower 2 bits
                //                ctmp+=mapAB2uc[aLine];
                if((j+1)%4==0) // full, plug in and refresh
                {
                    aLine.push_back(ctmp);
                    ctmp=0;
                }
            }
            if(ctmp!=0) // somthing-left non-full
            {
                // move them to the left
                ctmp <<= 2*(4-(j%4));
                aLine.push_back(ctmp);
            }
                //             translate               base         + offset
//                aLine = aLine + "\t" + vecUC2STR[ vecGenotypes[j*numberSNPs + i     ] ];
            mapSNP2Line[stmp]=aLine;
        }
        progClear();

        vector<unsigned char>().swap(vecGenotypes);
        map<string, string>().swap(mapSNPName2Pos);
        map<string, string>().swap(mapSNPName2Chr);
        // sort SNP
        sort(vecSNPinMap.begin(), vecSNPinMap.end(), SNP::SNPcomp);

        // **************** start writing
        L("s", 1, "Start writing");
        string outputFileName(abgFileName+".vcf");
        ofstream ofs (outputFileName.c_str());
        if(!ofs)
        {
            E("sss", "Initializing output file {", outputFileName.c_str(), "} failed. Aborting. X");
            exit(-10);
        }

        L("s", 2, "Writing out the header");
        ofs <<
        "##fileformat=VCFv4.1" << endl <<
        "##INFO=<ID=BeagleInputFile,Description=\"Generated by SNPipeline.pooledAB2VCF41\">" << endl <<
        "##INFO=<ID=DESC,Description=\"REF and ALT were given arbitrary from AA, AB and BB coding system\">" << endl <<
        "##REF=<ID=G,Description=\"Arbitrary set\">" << endl <<
        "##ALT=<ID=A,Description=\"Arbitrary set\">" << endl <<
        "##INFO=<ID=Genotype,Number=-10,Type=Integer,Description=\"AA, 0/0\">" << endl <<
        "##INFO=<ID=Genotype,Number=0,Type=Integer,Description=\"AB, 0/1\">" << endl <<
        "##INFO=<ID=Genotype,Number=10,Type=Integer,Description=\"BB, 1/1\">" << endl <<
        "##INFO=<ID=Genotype,Number=6.5,Type=Float,Description=\"--, ./.\">" << endl <<
        "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">" << endl <<
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for(i=0; i<numberSamples; ++i)
            ofs << "\t" << vecSampleIDs[i];
        ofs << endl;

        L("s", 2, "Writing out the genotypes for each SNP");
        // foreach SNP as an output line
        tn = mapSNP2Line.size()*.01;
        if(tn<100)
            tn=100;
        const string RAQPIF ("\tG\tA\t.\tPASS\t.\tGT");          // REF ALT QUAL. PASS INFO. FORMAT
        // vecSNPinMap is in order
        unsigned long prevPos(vecSNPinMap[0].pos+100);  // make 1st value unequal
        for(i=0; i<vecSNPinMap.size(); ++i)
        {
            if(i%tn==0)
                progShow('W', i/tn);

            stmp=vecSNPinMap[i].ID;  // vecSNPNames[0]=="SNP_IDs"
            if(mapSNP2Line.find(stmp)==mapSNP2Line.end())   // not found, skip
                continue;

            if (vecSNPinMap[i].pos==prevPos)   // same position, duplicated SNP, skip
                continue;

            prevPos=vecSNPinMap[i].pos;

            // print meta info #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
            ofs << vecSNPinMap[i].chr        // CHROM
            << "\t" << vecSNPinMap[i].pos       // POS
            << "\t" << stmp                       // ID
            << RAQPIF;         // REF ALT QUAL. PASS INFO. FORMAT
            // foreach animal as genotype output column
            aLine=mapSNP2Line[stmp];
            for(j=0; j<numberSamples; ++j)
            {
                // j   is genotype index
                // j/4 is char index
                gtSum = 2 * (3 - (j % 4));        // shift
                ofs << "\t" << vecUC2STR[(( (3 << gtSum) & aLine[j/4] ) >> gtSum)];
            }
            ofs << endl;
        }
        progClear();
        ofs.close();

        if(setDiscardedSNP.size()!=0)
        {
            outputFileName=abgFileName+".dscNM";
            ofstream ofs (outputFileName.c_str());
            I("iss", 1, setDiscardedSNP.size(), "SNPs not in map file were discarded. Detailed in file", outputFileName.c_str());
            for(set<string>::iterator sit=setDiscardedSNP.begin(); sit!=setDiscardedSNP.end(); ++sit)
                ofs << (*sit) << endl;
            ofs.close();
        }

        if(setDiscardedSNPMS.size()!=0)
        {
            outputFileName=abgFileName+".dscMS";
            ofstream ofs (outputFileName.c_str());
            I("iss", 1, setDiscardedSNPMS.size(), "no call SNPs were discarded. Detailed in file", outputFileName.c_str());
            for(set<string>::iterator sit=setDiscardedSNPMS.begin(); sit!=setDiscardedSNPMS.end(); ++sit)
                ofs << (*sit) << endl;
            ofs.close();
        }
        L("s", 1, "Done.");
    }
}

int main(const int argc, const char * argv[])
{
    libcbk::SNPipeline_pooledAB2VCF41(argc, argv);
    return 0;
}
