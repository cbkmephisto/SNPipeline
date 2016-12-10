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
//  SNPipeline.abgExcl.cpp
//  libHailins
//
//  Created by Hailin SU on 3/3/14.
//

#include<iostream>
#include<fstream>
#include<sstream>
#include<set>
#include<vector>
#include<cstring>

using namespace std;

int main (int argc, char * const argv[])
{
	const string version ("    SNPipeline.abgExcl    Version<0.01>    Mar 2014    hailins@iastate.edu");
    const string desc    ("        Pipe the ab-genotype file to a new one excluding animalIDs in the blackLisk file.");
	cout << endl << version << endl << endl << desc << endl << endl;

    if (argc != 5 || strlen(argv[3])!=2 || (argv[3][0] != '-' && argv[3][1] != 'o'))
    {
		cout << "Syntax: "<<argv[0]<<" <blackList_file_name> <ab-genotype_file_name> -o <abg-output_file_name> " << endl;
		return(1);
	}

    const string blackList(argv[1]), inputFileName(argv[2]), outputFileName(argv[4]);
    ifstream infs (blackList.c_str());                              //check for and open the ascii file genotype file
    if(!infs)
    {
        cerr << " # Couldn't open blackList file: " << blackList.c_str() << endl;
        return (-1);
    }
    set<string> setBlackList;
    string stmp;

    // loading blackList
    cout << " > Loading blackList ..." << endl;
    while(infs >> stmp)
        setBlackList.insert(stmp);
    infs.close();
    cout << "   - " << setBlackList.size() << " animal IDs in the blackList" << endl;
    infs.open(inputFileName.c_str());

    if(!infs)
    {
        cerr << " # Couldn't open ab-genotype data file: " << inputFileName.c_str() << endl;
        return (-1);
    }

    ofstream ofs;

    string aLine;

    // go header line
    stringstream sst;
    unsigned rec(0), recb(0), reco(0);

    if(getline(infs, aLine))
    {
        if(aLine[0]!='#' || aLine[8]!='#' || aLine[4]!='_')
        {
            cerr << " # NOT A VALID AB-GENOTYPE FILE." << endl << endl;
            return -1;
        }
        ofs.open(outputFileName.c_str());            //check for and open the binary genotype file
        if(!ofs)
        {
            infs.close();
            cerr << " # Couldn't open binary output file: " << outputFileName.c_str() << endl;
            return (-2);
        }

        ofs << aLine << endl;   // header line

        // data lines
        cout << " > Processing ..." << endl;
        while (getline(infs,aLine))
        {
            if(++rec%100==0)
            {
                cout << "\t" << rec << flush;
                if(rec%1500==0)
                    cout << endl;
            }
            sst.clear();
            sst << aLine;
            sst >> stmp;
            sst.str("");
            if(setBlackList.find(stmp)==setBlackList.end()) // write out not bl
            {
                ofs << aLine << endl;
                ++reco;
            }
            else
                ++recb;
        }
        cout << endl
        << " > Processed " << rec << " records." <<endl
        << "   - " << recb << " records excluded." <<endl
        << "   - " << reco << " records writen into [" << outputFileName << "]." <<endl;

    }
    return 0;
}