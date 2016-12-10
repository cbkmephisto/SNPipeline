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
//  SNPipeline.abg2bin.cpp
//  libHailins
//
//  Created by Hailin SU on 2/27/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "Eigen/Dense"
#include "../splitablestring.h"

using namespace std;
using namespace Eigen;

int main (int argc, char * const argv[])
{
    const string version ("    abg2bin    Version<0.02>    Apr 2015    hailins@iastate.edu");
    
    // Version<0.01>    Feb 2014
    // Version<0.02>    Apr 2015: adding a2g and g2b options
    cout << endl << version << endl << endl;
    
    
    string inputFileName, outputFileName;
    //    if (argc == 1)
    //    {
    //        cout << " > Enter genotype file name\n";
    //        cin >> inputFileName;
    //    }
    //    else if (argc == 2)
    if(argc == 2)
    {
        cout << " > Opening "<< argv[1]<<" as ab-genotype file to convert to binary"<<endl;
        inputFileName = argv[1];
    }
    else if (argc == 3 && (string(argv[1])=="-a2g" || string(argv[1])=="-g2b"))
    {
        cout << " > Opening "<< argv[2]<<" as input file to convert: " << argv[1] <<endl;
        inputFileName = argv[2];
    }
    else
    {
        cout << "Syntax: "<<argv[0]<<" [-a2g | -g2b] genotype_file_name " << endl;
        cout << "   * option -a2g converts the ab-genotype file to a -10/0/10 coded genotype file, replace missing data with column mean" << endl;
        cout << "      - non-segregating SNP locus will be coded as -77" << endl;
        cout << "   * option -g2b converts the -10/0/10 genotype file to a binary file" << endl;
        cout << "   * if neither -a2g nor -g2b was specified, abg2bin will convert ab-genotype file to binary format coding missing value as \"1\"" << endl;
        return (-1);
    }
    
    /*********************************************************************
     *********************************************************************
     *  abg2bin
     *********************************************************************
     *********************************************************************/
    if (argc == 2)  // abg2bin, original
    {
        outputFileName = inputFileName + ".newbin";
        
        ifstream infs (inputFileName.c_str());                              //check for and open the ascii file genotype file
        
        if(!infs)
        {
            cerr << " # Couldn't open ab-genotype data file: " << inputFileName.c_str() << endl;
            return (-1);
        }
        
        ofstream ofs(outputFileName.c_str(), ofstream::binary);            //check for and open the binary genotype file
        if(!ofs)
        {
            infs.close();
            cerr << " # Couldn't open binary output file: " << outputFileName.c_str() << endl;
            return (-1);
        }
        
        splitablestring aLine;
        
        // go header line
        vector<string> sx;
        unsigned numberMarkers, i, rec(0);
        map<string, short> umgc;
        
        if(getline(infs, aLine))
        {
            //  0123456789
            // "#SNP_IDs#ID1 ID2 ID3 ..."
            for(i=0; i<9; ++i)
                aLine[i]=' ';
            //  0123456789
            // "         ID1#ID2#ID3#..."
            for(i=9; i<aLine.length(); ++i)
                if(aLine[i]=='#')
                    aLine[i]=' ';
            
            sx=aLine.split(' ');
            numberMarkers=unsigned(sx.size());
            cout << " +--> " << numberMarkers << " markers detected." << endl;
            //        cout << aLine << endl;
            string headerLine ("SNP_IDs.bt" + aLine);
            
            umgc["AA"]=-10;
            umgc["AB"]=0;
            umgc["BB"]=10;
            umgc["--"]=1;
            umgc["ST"]=1;
            
            short rowi[numberMarkers];      // this is the version of the marker covariates stored as real numbers
            
            // write out the number of markers as an unsigned integer
            ofs.write((char*) &numberMarkers, sizeof(numberMarkers));    // once we know how many bytes in an unsigned integer
            
            unsigned hdrlength ((unsigned)headerLine.length());
            ofs.write((char*) &hdrlength, sizeof(hdrlength));        // write out the length of the header
            ofs.write((char*) headerLine.c_str(), hdrlength);                // write out the character header of SNP names itself
            
            const unsigned markerLength (numberMarkers*sizeof(short));
            
            while (getline(infs,aLine))
            {
                if(++rec%100==0)
                {
                    cout << "\t" << rec << flush;
                    if(rec%700==0)
                        cout << endl;
                }
                
                sx=aLine.split(' ');
                
                // first column  (column 0) is animal id
                hdrlength = unsigned(sx[0].length());
                ofs.write((char*) &hdrlength, sizeof(hdrlength));       // write out the length of the animalID
                
                ofs.write((char*) sx[0].c_str(), hdrlength);        // write out the animalID itself
                
                for(i=0; i<numberMarkers; ++i)
                    rowi[i] = umgc[sx[i+1]];
                
                ofs.write((char*) rowi, markerLength);		  // write out the actual marker genotypes
            }
            cout << endl << " > Processed " << rec << " records." <<endl;
        }
        return 0;
    }
    /*********************************************************************
     *********************************************************************
     * abg2gen, 'AA'->-10, 'BB'->10, 'AB'->0, '--'->mean
     *********************************************************************
     *********************************************************************/
    else if (argv[1][3]=='g') // abg2gen
    {
        outputFileName = inputFileName + ".abg2gen";
        string ofn2 = inputFileName + ".colmean";
        
        ifstream infs (inputFileName.c_str());                              //check for and open the ascii file genotype file
        
        if(!infs)
        {
            cerr << " # Couldn't open ab-genotype data file: " << inputFileName.c_str() << endl;
            return (-1);
        }
        
        ofstream ofs(outputFileName.c_str());            //check for and open the binary genotype file
        if(!ofs)
        {
            infs.close();
            cerr << " # Couldn't open output file: " << outputFileName.c_str() << endl;
            return (-1);
        }
        
        ofstream ofscm(ofn2.c_str());
        
        splitablestring aLine;
        
        // go header line
        vector<string> sx;
        unsigned numberMarkers(0), numberAnimals(0), i, j;
        map<string, short> umgc;
        
        // read the 1st time to get numberAnimals
        if(getline(infs, aLine))    // the headerLine
        {
            while(getline(infs, aLine))
                ++numberAnimals;
        }
        
        infs.close();
        infs.open(inputFileName.c_str());
        // read the 2nd time to read in the genotype matrix
        if(getline(infs, aLine))
        {
            //  0123456789
            // "#SNP_IDs#ID1 ID2 ID3 ..."
            string headerLine (aLine.c_str());
            
            vector<unsigned> vecN;
            vector<float> vecSum;
            
            for(i=0; i<9; ++i)
                aLine[i]=' ';
            //  0123456789
            // "         ID1#ID2#ID3#..."
            for(i=9; i<aLine.length(); ++i)
                if(aLine[i]=='#')
                    aLine[i]=' ';
            
            sx=aLine.split(' ');
            numberMarkers=unsigned(sx.size());
            vecN.resize(numberMarkers);
            vecSum.resize(numberMarkers);
            for(i=0; i<numberMarkers; ++i)
            {
                vecN[i]=numberAnimals;
                vecSum[i]=0.0;
            }
            cout << " +--> " << numberMarkers << " markers detected." << endl;
            //        cout << aLine << endl;
            vector<string>().swap(sx);
            
            umgc["AA"]=-10;
            umgc["AB"]=0;
            umgc["BB"]=10;
            umgc["--"]=99;   // trick: put missing as 99, and subtract 1 from numberAnimals for that column
            umgc["ST"]=99;   // trick: put missing as 99, and subtract 1 from numberAnimals for that column
            
            MatrixXi    X;                  //  the incident matrix
            X.resize(numberAnimals, numberMarkers);
            X.setZero();
            unsigned nMissing(0);
            vector<string> vecAnmIDs;
            for(i=0; i<numberAnimals; ++i)
            {
                if(i>0 && i%100==0)
                {
                    cout << "\t" << i << flush;
                    if(i%700==0)
                        cout << endl;
                }
                
                getline(infs, aLine);
                sx=aLine.split();
                
                sx=aLine.split(' ');
                // first column  (column 0) is animal id
                vecAnmIDs.push_back(sx[0]);
                
                for(j=0; j<numberMarkers; ++j)
                {
                    X(i, j) = umgc[sx[j+1]];
                    if(sx[j+1]=="--")           // missing
                        --vecN[j];
                    else
                        vecSum[j]+=X(i, j);
                }
            }
            cout << endl << " > Loaded " << numberAnimals << " records." <<endl;
            cout << " > Filling missing genotypes..." <<endl;
            ofs << headerLine << endl;
            // calc the mean and fill missing in X
            for(j=0; j<numberMarkers; ++j)
            {
                if(vecN[j]>0)
                    vecSum[j]/=vecN[j];
                else
                    vecSum[j]=-77;        // no mean: genotype missing the whole column
                
                if(j!=0)
                    ofscm << " " << vecSum[j];
                else
                    ofscm << vecSum[j];
            }
            ofscm << endl;
            // now vecSum[j] is the mean for col j
            for(i=0; i<numberAnimals; ++i)
            {
                if(i>0 && i%100==0)
                {
                    cout << "\t" << i << flush;
                    if(i%700==0)
                        cout << endl;
                }
                
                ofs << vecAnmIDs[i];
                for(j=0; j<numberMarkers; ++j)
                {
                    if(X(i, j)==99)
                    {
                        X(i, j)=vecSum[j];
                        ++nMissing;
                    }
                    ofs << " " << X(i, j);
                }
                ofs << endl;
            }
            ofs.close();
            ofscm.close();
            cout << endl << " > " << nMissing <<" (" << int(1000*nMissing/(numberAnimals*numberMarkers))/10.0<< "%) missing genotypes filled by column mean." <<endl;
        }
        return 0;
    }
    /*********************************************************************
     *********************************************************************
     * gen2bin
     *********************************************************************
     *********************************************************************/
    else                       // gen2bin
    {
        outputFileName = inputFileName + ".newbin";
        
        ifstream infs (inputFileName.c_str());                              //check for and open the ascii file genotype file
        
        if(!infs)
        {
            cerr << " # Couldn't open -10/0/10-genotype data file: " << inputFileName.c_str() << endl;
            return (-1);
        }
        
        ofstream ofs(outputFileName.c_str(), ofstream::binary);            //check for and open the binary genotype file
        if(!ofs)
        {
            infs.close();
            cerr << " # Couldn't open binary output file: " << outputFileName.c_str() << endl;
            return (-1);
        }
        
        splitablestring aLine;
        
        // go header line
        vector<string> sx;
        unsigned numberMarkers, i, rec(0);
        map<string, short> umgc;
        
        if(getline(infs, aLine))
        {
            //  0123456789
            // "#SNP_IDs#ID1 ID2 ID3 ..."
            for(i=0; i<9; ++i)
                aLine[i]=' ';
            //  0123456789
            // "         ID1#ID2#ID3#..."
            for(i=9; i<aLine.length(); ++i)
                if(aLine[i]=='#')
                    aLine[i]=' ';
            
            sx=aLine.split(' ');
            numberMarkers=unsigned(sx.size());
            cout << " +--> " << numberMarkers << " markers detected." << endl;
            //        cout << aLine << endl;
            string headerLine ("SNP_IDs.bt" + aLine);
            
            short rowi[numberMarkers], genotype;      // this is the version of the marker covariates stored as real numbers
            stringstream sst;
            
            // write out the number of markers as an unsigned integer
            ofs.write((char*) &numberMarkers, sizeof(numberMarkers));    // once we know how many bytes in an unsigned integer
            
            unsigned hdrlength ((unsigned)headerLine.length());
            ofs.write((char*) &hdrlength, sizeof(hdrlength));        // write out the length of the header
            ofs.write((char*) headerLine.c_str(), hdrlength);                // write out the character header of SNP names itself
            
            const unsigned markerLength (numberMarkers*sizeof(short));
            
            while (getline(infs,aLine))
            {
                if(++rec%100==0)
                {
                    cout << "\t" << rec << flush;
                    if(rec%700==0)
                        cout << endl;
                }
                
                sx=aLine.split(' ');
                
                // first column  (column 0) is animal id
                hdrlength = unsigned(sx[0].length());
                ofs.write((char*) &hdrlength, sizeof(hdrlength));       // write out the length of the animalID
                
                ofs.write((char*) sx[0].c_str(), hdrlength);        // write out the animalID itself
                
                for(i=0; i<numberMarkers; ++i)
                {
                    sst << sx[i+1];
                    sst >> genotype;
                    sst.clear();
                    rowi[i] = genotype;
                }
                ofs.write((char*) rowi, markerLength);		  // write out the actual marker genotypes
            }
            cout << endl << " > Processed " << rec << " records." <<endl;
        }
        return 0;
    }
}