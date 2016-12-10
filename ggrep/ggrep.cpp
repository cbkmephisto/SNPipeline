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
//  ggrep.cpp
//  libHailins
//
//  Created by Hailin SU on 3/26/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#include "ggrep.h"

using namespace std;

void usage(char* argv0)
{
    cout << endl;
    cout << "ggrep: Group Grep by hailins@iastate.edu Wed Mar 26 2:18PM 2014" << endl << endl;
    cout << "Usage:" << endl << endl
    << argv0 << " lst_file target_file > newFile" << endl
    << argv0 << " -v blacklst_file target_file  > newFile" << endl << endl;
    cout << "This program" << endl
    << " - greps the lines acording to that 1st col of target_file listed in lst_file, OR" << endl
    << " - excludes the lines of target_file for all the 1st col listed in blacklst_file (-v)." << endl
    << " - lines started with # will always be grepped and printed out." << endl
    << endl;
}

int main (const int argc, char** argv)
{
    if( (argc!=3 && argc!=4) || (argc==4 && string(argv[1])!="-v") )
    {
        usage(argv[0]);
        exit (-1);
    }

    bool v(false);
    string lstFile, tgtFile;
    if(argc==3)
    {
        lstFile=argv[1];
        tgtFile=argv[2];
    }
    else
    {
        v=true;
        lstFile=argv[2];
        tgtFile=argv[3];
    }

    // read the list
    ifstream infs;
    infs.open(lstFile.c_str());
    if(!infs)
    {
        cerr << "List file " << lstFile.c_str() << " could not be open to read." << endl;
        exit (-2);
    }

    set<string> setRef;
    string stmp, aLine;
    while(infs >> stmp)
        setRef.insert(stmp);

    infs.close();

    // pipe the target
    infs.open(tgtFile.c_str());
    if(!infs)
    {
        cerr << "Target file " << tgtFile.c_str() << " could not be open to read." << endl;
        exit (-3);
    }

    bool inList;
    stringstream sst;
    while(getline(infs, aLine))
    {
        if(aLine[0]=='#')   // reserve the header line for trimDown2Map
            inList=true;
        else
            inList=false;
        sst << aLine;
        sst >> stmp;
        sst.str("");
        sst.clear();

        if(setRef.find(stmp)!=setRef.end()) // in the list
        {
            inList=true;
//            setRef.erase(stmp);
        }
        if((!v && inList) || (v && !inList))
           cout << aLine << endl;
    }
    infs.close();
}
