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
//  abg2M.cpp
//  libHailins
//
//  Created by Hailin Su on 12/14/15.
//

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "../splitablestring.h"

using namespace std;

int main (int argc, char * const argv[])
{
    const string version ("    abg2M    Version<0.01>    Dec 2015    hailins@iastate.edu");
    
    // Version<0.01>    Dec 2015
    cout << endl << version << endl << endl;
    

    if(argc == 4)
    {
        string  ruleFileName(argv[1]),
                abgFileName( argv[2]),
                outputPrefix(argv[3]),
                outputAID,
                outputSID,
                outputMTX;

        outputAID = outputPrefix+".anm_ID";
        outputSID = outputPrefix+".SNP_ID",
        outputMTX = outputPrefix+".mtx_ds";
        splitablestring aLine;
        vector<string>  vecTmp;
        string          vecCoding[4];

        unsigned utmp;

        map<string, string> mapCodingXref;

        ////////////////////////////////////////////////////////////////////////////////////////////
        // read rule
        ////////////////////////////////////////////////////////////////////////////////////////////
        cout << " > Opening "<< argv[1]<<" as rule file"<<endl;
        ifstream infs (ruleFileName.c_str());
        if(!infs)
        {
            cerr << " # Couldn't open rule file: " << ruleFileName.c_str() << endl;
            return (-1);
        }

        while(getline(infs, aLine))
        {
            vecTmp=aLine.split();
            if((utmp=unsigned(vecTmp.size()))<2)
                continue;
            mapCodingXref[vecTmp[0]]=vecTmp[utmp-1];
        }
        infs.close();
        
        if(mapCodingXref.find("AA")==mapCodingXref.end()
           || mapCodingXref.find("AB")==mapCodingXref.end()
           || mapCodingXref.find("BB")==mapCodingXref.end()
           || mapCodingXref.find("--")==mapCodingXref.end())
        {
            cerr << " # rule file did not specified rule(s) for AA, AB, BB, or --" << endl;
            return (-2);
        }

        // missing is now allowed!
        vecCoding[0]=mapCodingXref["AA"];
        vecCoding[1]=mapCodingXref["AB"];
        vecCoding[2]=mapCodingXref["BB"];
        vecCoding[3]=mapCodingXref["--"];

        ////////////////////////////////////////////////////////////////////////////////////////////
        // convert
        ////////////////////////////////////////////////////////////////////////////////////////////
        cout << " > Opening "<< argv[2]<<" as ab-genotype file to convert to marker matrix (dense!)"<<endl;

        infs.open(abgFileName.c_str());
        if(!infs)
        {
            cerr << " # Couldn't open ab-genotype data file: " << abgFileName.c_str() << endl;
            return (-1);
        }
        
        ofstream ofsAID(outputAID.c_str());
        ofstream ofsSID(outputSID.c_str());
        ofstream ofsMTX(outputMTX.c_str());
        if(!ofsAID || !ofsSID || !ofsMTX)
        {
            infs.close();
            cerr << " # Couldn't open files to write." << endl;
            return (-1);
        }
        
        // go header line, write out SID
        unsigned i, rec(0);
        unsigned short gtmp;
        bool invalidGT = false;
        string gtLine;
        if(getline(infs, aLine))
        {
            vecTmp=aLine.split('#');
            utmp=unsigned(vecTmp.size());
            for(i=1; i<utmp; ++i)
                ofsSID << vecTmp[i] << endl;
            
            while (getline(infs,aLine))
            {
                if(++rec%100==0)
                {
                    cout << "\t" << rec << flush;
                    if(rec%700==0)
                        cout << endl;
                }
                if(aLine[0]=='#')
                    continue;

                gtLine = "";
                vecTmp=aLine.split();
                // first column  (column 0) is animal id
                ofsAID << vecTmp[0] << endl;
                utmp=unsigned(vecTmp.size());
                for(i=1; i<utmp; ++i)
                {
                    gtmp=vecTmp[i][0]-'A'+vecTmp[i][1]-'A';
                    if (vecTmp[i][0] == vecTmp[i][1] && vecTmp[i][1] == '-')
                    {
                        gtmp = 3;
                        // cout << "gtmp is set to 3" << endl;
                    }
                    else if(gtmp>3)
                    {
                        cerr << " # Something wroing was in the ab-genotype coding: " << vecTmp[i] << " for animal "
                        << vecTmp[0] << endl;
                        invalidGT = true;
                        // cout << gtmp << " " << vecTmp[i][0] << " " << vecTmp[i][1] << (vecTmp[i][0] == vecTmp[i][1]) << ('-' == vecTmp[i][1]) << endl;
                        break;
                    }
                    if(gtLine.length())
                        gtLine += " ";
                    gtLine += vecCoding[gtmp];
//                    ofsMTX << vecCoding[gtmp] << " ";
                }
                ofsMTX << gtLine << endl;
                if(invalidGT)
                    break;
            }
            cout << endl << " > Processed " << rec << " records." <<endl;
        }
        infs.close();
        ofsSID.close();
        ofsAID.close();
        ofsMTX.close();
        return 0;
    }
    else
    {
        cout << "Syntax: "<<argv[0]<<" <rule_file> <ab-genotype_file_name> <output_prefix>" << endl;
        cout << " - <rule_file> is the file defining conversion. For a -10, 0, 10 coding system, use sample below:" << endl;
        cout << "**** file content below" << endl;
        cout << "AA -10" << endl;
        cout << "AB 0" << endl;
        cout << "BB 10" << endl;
        cout << "-- ??" << endl;
        cout << "**** file content above" << endl;
        cout << " - <output_prefix> defines the prefix name of the output files. 3 output files will be generated:" << endl;
        cout << "   - <output_prefix>.anm_ID, storing animal IDs, one per line" << endl;
        cout << "   - <output_prefix>.SNP_ID, storing SNP IDs, one per line" << endl;
        cout << "   - <output_prefix>.mtx_ds, storing the dense coded matrix" << endl;
        return (-1);
    }
}