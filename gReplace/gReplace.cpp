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
//  gReplace.cpp
//  libHailins
//
//  Created by Hailin SU on 4/10/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#include "gReplace.h"
#include "../splitablestring.h"

using namespace std;

void usage(char* argv0)
{
    cout << endl;
    cout << "gReplace: Group Replace by hailins@iastate.edu Thr Apr 10 10:58 AM 2014" << endl << endl;
    cout << "Usage:" << endl << endl
    << argv0 << " 2col_xref_file target_file > newFile" << endl << endl;
    cout << "This program" << endl
    << " - replaces the words in the col1 of 2col_xref_file to col2 if the col1 words appeared in the target_file" << endl
    << " - print out on screen, delimetered by whitespace." << endl
    << endl;
}

int main (const int argc, char** argv)
{
    if(argc!=3)
    {
        usage(argv[0]);
        exit (-1);
    }
    string lstFile(argv[1]), tgtFile(argv[2]);

    // read the list
    ifstream infs;
    infs.open(lstFile.c_str());
    if(!infs)
    {
        cerr << "Map file " << lstFile.c_str() << " could not be open to read." << endl;
        exit (-2);
    }

    map<string, string> mapXref;
    string stmp;
    splitablestring aLine;
    vector<string> vtmp;
    cerr << " - loading xref" << endl;
    while(getline(infs,aLine))
    {
        vtmp=aLine.split();
        if(vtmp.size()!=2)
        {
            cerr << "xref file " << lstFile.c_str() << " is not a 2-col-through file." << endl;
            cerr << "line content [" << aLine << "] contains " << vtmp.size() << " fields. aborting" << endl;
            infs.close();
            exit(-3);
        }
        mapXref[vtmp[0]]=vtmp[1];
    }
    infs.close();
    vector<string>().swap(vtmp);
    // pipe the target
    infs.open(tgtFile.c_str());
    if(!infs)
    {
        cerr << "Target file " << tgtFile.c_str() << " could not be open to read." << endl;
        exit (-3);
    }

//    for(map<string, string>::iterator mit=mapXref.begin(); mit!=mapXref.end(); ++mit)
//        cerr << mit->first << "   ->    " << mit->second << endl;
    
    cerr << " - x-referencing" << endl;
    stringstream sst;
    while(getline(infs, aLine))
    {
        sst << aLine;
        sst >> aLine;
        if(mapXref.find(aLine)!=mapXref.end())
            aLine=mapXref[aLine];

        while(sst>>stmp)
        {
            aLine+=" ";
            if(mapXref.find(stmp)!=mapXref.end()) // in the list
                aLine+=mapXref[stmp];
            else
                aLine+=stmp;
        }
        sst.str("");
        sst.clear();
        cout << aLine << endl;
    }
    infs.close();
}
