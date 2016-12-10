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
//  bout2genotype.cpp
//  libHailins
//
//  Created by Hailin SU on 11/4/13.
//  Copyright (c) 2013 iastate. All rights reserved.
//
//  bout2genotype

#include "../appInfo.h"
#include "geneticMap.h"
#include "beagleOutput.h"
#include <iomanip>

using namespace std;
namespace libcbk
{
    class bout2genotype : public appInfo
    {
    public:
        bout2genotype(int argc, const char * argv[]);
    private:
        void theMain(int argc, const char * argv[]);
    };

    bout2genotype::bout2genotype(int argc, const char * argv[])
    {
        progName    = "bout2genotype";
        version     = "2013-11-11 09:27 CDT";
        created     = "2013-11-04";
        updated     = "2014-06-30";
        descrip     = "read beagle output files to generate files for doing SNPRelate PCA.\n\tmatrix: 0-1-2 & SNP-in-row genotype files";
        displayInfo();
        lg.initLogger();
        theMain(argc, argv);
        lg.finalizeLogger();
        /*
         updates
         2013-11-04 initialized
         2013-11-06 output to 4 files: snp_id, animal_id, matrix, map
         2013-11-11 output in map position order; output every 10 SNPs
         */
    }

    void bout2genotype::theMain(int argc, const char * argv[])
    {
        makeSure(argc >= 4, argv[0], "incFile mapFile beagle_(*.phased)");

        justAmap jam;
        I("s", 1, "constructing genetic map via incFile and mapFile...");
        map<short, geneticMap> gmaps = geneticMap::constructGeneticMap(argv[1], argv[2], jam);

        ofstream ofsm ("bo2g.matrix"), ofsa("bo2g.aid"), ofss("bo2g.sid"), ofsp("bo2g.map");
        if(!ofsm || !ofsa || !ofss || !ofsp)
        {
            ofsm.close(); ofss.close(); ofsa.close(); ofsp.close();
            cout << "Do not have the write access in current folder. Terminated." << endl;
            exit (-150);
        }

//      map<string, string>::iterator mit; // for output-position order
        map<unsigned long, string>::iterator mit;
        unsigned genotype, allele;
        beagleOutput bout;
//        long every10=-1;
        for(unsigned i=3; i<argc; ++i)
        {
            I("sss", 1, "loading beagle output file : [", argv[i], "]");
            bout.load(argv[i], jam);
            if(i==3)    // 1st file, write out the anima IDs
            {
//                ofs << left << setw(25) << setfill(' ') << "SNP___ANM";
                for(unsigned j=0; j<bout.animalIDs.size(); ++j)
                    ofsa << bout.animalIDs[j] << "\t" << bout.breed << endl;
            }
            stringstream ss;                // declaring outside the loop may cause insufficient memory

            I("s", 1, "writing");
/*            // for each SNP
            for(mit=bout.mapMKN2Line.begin(); mit!=bout.mapMKN2Line.end(); ++mit)
            {
                //  mit->first = SNPid, mit->second = SNP genotype Line
//                ofs << left << setw(25) << setfill(' ') << mit->first << " ";       // 1st column: SNP id
                ofss << mit->first << endl;       // 1st column: SNP id
                ofsp << jam.mapMKN2line[mit->first] << endl;
                ss << mit->second;  // load
                while(ss>>genotype)
                {
                    ss>>allele;
                    genotype+=allele;
                    ofsm << genotype << " ";       // columns: genotype
                }
                ofsm << endl;
                ss.clear();
            }
 */// 2013-11-11 output in map position order + output every 10 SNP
            // mit = iterator of     map<unsigned long, string>  mapPos2Name;  //  pos -> markerName
            // mit->first = position, mit->second = SNP_id       in the map
//      2014-06-30 filter out unsegragating snps
            string genLine (""), genChar;
            bool bSeg(false);
            unsigned lastG(255);
            stringstream ssg;
            for(mit=gmaps[bout.chr].mapPos2Name.begin(); mit!=gmaps[bout.chr].mapPos2Name.end(); ++mit)
            {
                // see if the current SNP is in the list: because of the include file
                if(bout.mapMKN2Line.find(mit->second)==bout.mapMKN2Line.end())  // not found
                    continue;
                // get every 10 SNP
//                if(++every10 % 10 != 0) continue;

                ss << bout.mapMKN2Line[mit->second];  // load the genotype
                while(ss>>genotype)
                {
                    ss>>allele;
                    genotype+=allele;
                    if(lastG==255)
                        lastG=genotype;
                    else if (!bSeg && lastG!=genotype)
                        bSeg=true;  // set the flag
                    ssg << genotype << " ";
//                    ofsm << genotype << " ";       // columns: genotype
                }

                // only output segragating snps
                if(bSeg)
                {
                    while(ssg >> genChar)
                    {
                        genLine.append(genChar);
                        genLine.append(" ");
                    }
                    ofsm << genLine << endl;
                    ofss << mit->second << endl;                                    // write out the SNP id
                    ofsp << jam.mapMKN2line[mit->second] << endl; // write out the line for the current SNP of the map file
                }

                genLine="";
                lastG=255;
                bSeg=false;
                ss.clear();
                ssg.clear();
            }
            bout.clear();
        }
        ofsm.close();
        ofsa.close();
        ofss.close();
        ofsp.close();
    }
}

int main(int argc, const char * argv[])
{
    libcbk::bout2genotype(argc, argv);
    return 0;
}