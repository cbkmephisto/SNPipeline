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
//  beagleOutput.cpp
//  haplotype_diversity
//
//  Created by Hailin SU on 9/11/13.
//  Copyright (c) 2013 iastate. All rights reserved.
//
//  updated 2013-10-01
//      I noticed that I've made a big mistake: I should not use a mapcount2hl, because the reverse-mapping couldn't correctly map from
//          count -> hl: many hl's would have the same count
//
#include "beagleOutput.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void beagleOutput::load(string fileName, justAmap &jmp)
{
    ifstream ifb (fileName.c_str());
    if(!ifb)
    {
        cout << "XXXXXXXXXXXX bealge output file not found: " << fileName << endl;
        exit(-20);
    }
    /***********
                    filling out the header: animalIDs
     ***********/
    cout << "                + filling out the header: animalIDs..." << endl;
    clear();
    string aLine;
    stringstream sst;
    ifb >> aLine >> aLine;  // the 'I' and 'Animal'
    getline(ifb, aLine);
    sst<<aLine;
    // first 2 columns, get the breed
    if(sst>>aLine)
    {
        breed=aLine.substr(0,3);    // get the breed
        sst>>aLine;
        animalIDs.push_back(aLine);
    }
    while(sst>>aLine)       // the animalID appears twice
    {
        sst>>aLine;
        animalIDs.push_back(aLine);
    }
    sst.clear();
    haplotypeRC.resize(animalIDs.size()*2);     //denotes haplotype is rare/common
    for(unsigned i=0; i<haplotypeRC.size(); ++i)
        haplotypeRC[i]="";

    cout << "                   - " << animalIDs.size() << " reads" << endl;

    /***********
                    filling out the map: mapMKN2Line
     ***********/
    cout << "                + filling out the map   : mapMKN2Line..." << endl;
    string markerName;
    chr=0;
    while(ifb >> aLine)       // the 'M' flag
    {
        ifb >> markerName;
        if(jmp.mapMKN2chr.find(markerName)==jmp.mapMKN2chr.end())// not in the inc list
            continue;
        if(chr==0)
        {
            // first get the current chr from map
            chr = jmp.mapMKN2chr[markerName];
            cout << "                   - chr " << chr << endl;
        }
        getline(ifb, aLine);
        mapMKN2Line[markerName] = aLine;
    }
    ifb.close();
}

void beagleOutput::generateResult(geneticMap &gmp)
{
    generateResultW(gmp, (float)1.0);
}

void beagleOutput::generateResultW(geneticMap &gmp, float winSizeM)
{
    static string outFile1, outFile2, outFile3, outFile4;
    string strchr;
    stringstream sst;
    sst << chr;
    sst >> strchr;
    sst.clear();

    // "hdout.<breed>.<chr>.totl", counts number of haplotype alleles in each window
    outFile1 = "hdout." + breed + "." + strchr + ".totl";
    ofstream ofstotl (outFile1.c_str());

    // "hdout.<breed>.<chr>.dtld", counts times of occurrance of each haplotype allele in each window
    outFile2 = "hdout." + breed + "." + strchr + ".dtld";
    ofstream ofsdtld (outFile2.c_str());

    // "hdout.<breed>.<chr>.racm", shows the structure of (ra)re / (c)o(m)mon haplotypes for each animal: 2 lines per animal.
    outFile3 = "hdout." + breed + "." + strchr + ".racm";
    ofstream ofsracm (outFile3.c_str());

    // "hdout.<breed>.<chr>.rprp", rare proportion
    outFile4 = "hdout." + breed + "." + strchr + ".rprp";
    ofstream ofsrprp (outFile4.c_str());

    if(!ofstotl || !ofsdtld || !ofsracm || !ofsrprp)
    {
        ofstotl.close(); ofsdtld.close(); ofsracm.close(); ofsrprp.close();
        cout << "XXXXXXXXXXXX cannot open files to write into." << endl;
        exit (-30);
    }
    unsigned long       curPOS (0);
    unsigned int        curWIN (0);
    static const unsigned long WINDOW_SIZE (winSizeM*1000000L);  // 1M
    string              theLine, lastLine("");
    //  go through the map (within the chr)
    map<unsigned long, string>::iterator mit;
    vector<string> queue;

    ofstotl << left
            << setw(23) << "number_of_diff_haps"
            << setw(16) << "prop(noh/allh)"
            << setw(20) << "mean_of_hap_occs"
            << setw(8)  << "chr"
            << setw(8)  << "window"
            << setw(8)  << "breed"
            << endl;
    ofsdtld << left
            << setw(19) << "hap_length"
            << setw(winSizeM*38) << "haplotype "
            << setw(8)  << "occr"
            << setw(16)  << "freq"
            << setw(8)  << "chr"
            << setw(8)  << "window"
            << setw(8)  << "breed"
            << setw(10)  << "homo_occr"
            << setw(10)  << "homo_freq"
            << endl;

    map<unsigned long, string> *theMap (&(gmp.mapPos2Name));
    for(mit = theMap->begin(); mit != theMap->end(); ++mit)
    {
        curPOS  = mit->first;               // first  = pos in map
        theLine = mapMKN2Line[mit->second]; // second = markerName mapped with pos, theLine = the line of beagle output containing given markerName
        while((1+curWIN)*WINDOW_SIZE<curPOS)
        {
            if(queue.size()==0)
            {
                ofstotl << left
                        << setw(23) << 0
                        << setw(16) << 0
                        << setw(20) << 0
                        << setw(8)  << chr
                        << setw(8)  << curWIN
                        << setw(8)  << breed
                        << endl;
//                lastLine=""; // in case of 0 SNP window[s]
            }
            else
            {
                eligantCalcAndWrite(ofstotl, ofsdtld, queue, curWIN, winSizeM*38);
                // storing last SNP of previous window
                lastLine=queue[queue.size()-1];
            }
            vector<string>().swap(queue);
            // 20150603: fill the gap between two windows: haplotype starting from last SNP of previous window
            if(lastLine.length())
                queue.push_back(lastLine);
            curWIN++;
        }
        queue.push_back(theLine);
    }
    if(queue.size()!=0)
    {
        eligantCalcAndWrite(ofstotl, ofsdtld, queue, curWIN, winSizeM*38);
        vector<string>().swap(queue);
    }
    cout << "                + output file [" << outFile1.c_str() << "] generated." << endl;
    cout << "                + output file [" << outFile2.c_str() << "] generated." << endl;

    unsigned r_count(0), c_count(0), j;
    for(unsigned i=0; i<haplotypeRC.size(); ++i) // write out the haplotype rare/common info
    {
        ofsracm << left << setw(24) << animalIDs[i/2] << haplotypeRC[i] << endl;
        for(j=0; j<haplotypeRC[i].length(); ++j)
        {
            if(haplotypeRC[i][j] == 'r')
                ++r_count;
            else if(haplotypeRC[i][j] == 'c')
                ++c_count;
        }
        if(i && i%2==0)
        {
            ofsrprp << left << setw(8) << breed << left << setw(8) << strchr << setw(8) << r_count/float(r_count+c_count) << endl;
            c_count = r_count = 0;
        }
    }
    cout << "                + output file [" << outFile3.c_str() << "] generated." << endl;
    cout << "                + output file [" << outFile4.c_str() << "] generated." << endl;
    //    cout << endl;
    ofstotl.close();
    ofsdtld.close();
    ofsracm.close();
    ofsrprp.close();
//    animalIDs.clear();
    vector<string>().swap(animalIDs);
//    mapMKN2Line.clear();
    map<string, string>().swap(mapMKN2Line);
    vector<string>().swap(haplotypeRC);
}
/*
string beagleOutput::getAnimalIDbyIndex(unsigned int ind)
{
    return animalIDs[ind/2];
}
*/

void beagleOutput::eligantCalcAndWrite(ofstream &ofstotl, ofstream &ofsdtld, vector<string> &queue, unsigned int curWIN, unsigned hapWidth)
{
    static unsigned i, j, homo_count;
    static unsigned long nHAPs(0), nSNPs(0);
//    static unsigned short ctmp;
    static string stmp;
    static stringstream sst;
    static vector<string> hapsl;    // string supports no limit
    static map<string, unsigned short> maphl2count;
    static map<string, unsigned short>::iterator mit;
    static float meanx, haplotype_frequency, homo_freq;
    static map<string, string> maphl2rareFlag;
    static map<string, unsigned short> maphl2homo_count;
/*                      A1  A1  A2  A2  A3  A3 ...
 *  queue[0]    SNP1
 *  queue[1]    SNP2
 *  queue[2]    SNP3
 *  ...         ...
 */
//    cout << "."; cout.flush();
    //  form the haplotypes
    nSNPs = queue.size();   // numberSNPs in the current window
    sst.clear();
    sst << queue[0];
    while(sst >> stmp)
        hapsl.push_back(stmp);  // This loop push_back the 1st SNP to all the hapsl
/*                  S1  S2  S3  ...
 *  hapsl[0]    H1
 *  hapsl[1]    H2
 *  hapsl[2]    H3
 *  ...         ...
 */
    for(i=1; i<nSNPs; i++)      // And this loop append all other SNPs in the queue to the hapsl by shifting left and add
    {
        sst.clear();
        sst << queue[i];
        j=0;
        while(sst >> stmp)
        {
            hapsl[j++] += stmp;
//            hapsl[j] |= ctmp;
//            j++;
        }
    }
    nHAPs=hapsl.size();
    // counting...
    for(i=0; i<nHAPs; i++)
    {
        stmp=hapsl[i];
        if(maphl2count.find(stmp)==maphl2count.end())    // new haplotype
            maphl2count[stmp] = 1;
        else                                                 // existed
            maphl2count[stmp] += 1;
        if(i%2) // check if homo every odd index
            if(stmp==hapsl[i-1])
            {
                if(maphl2homo_count.find(stmp)==maphl2homo_count.end())     // new homo
                    maphl2homo_count[stmp] = 1;
                else                                                            // existed
                    maphl2homo_count[stmp] += 1;
            }
    }
    meanx=0.0f;
    for (mit=maphl2count.begin(); mit!=maphl2count.end(); ++mit)
        meanx+=mit->second;
    meanx/=(float)nSNPs;
    ofstotl << left
            << setw(23) << maphl2count.size()
            << setw(16) << maphl2count.size() / (float)nHAPs
            << setw(20) << meanx
            << setw(8)  << chr
            << setw(8)  << curWIN
            << setw(8)  << breed
            << endl;

    for (mit=maphl2count.begin(); mit!=maphl2count.end(); ++mit)
    {
        ofsdtld.setf(ios::uppercase);
        ofsdtld << left << setw(19) << setfill(' ') << nSNPs;
        ofsdtld.unsetf(ios::uppercase);

        stmp=mit->first;    // the hap string
        if(maphl2homo_count.find(stmp)==maphl2homo_count.end())                // not in the homo map
            homo_count = 0;
        else
            homo_count = maphl2homo_count[stmp];
        homo_freq = 2.0 * homo_count / nHAPs;

        haplotype_frequency = mit->second/(float)nHAPs;

        ofsdtld << setw(hapWidth) << stmp << " "
                << setw(8)  << mit->second
                << setw(16) << haplotype_frequency
                << setw(8)  << chr
                << setw(8)  << curWIN
                << setw(8)  << breed
                << setw(10) << homo_count
                << setw(10) << homo_freq
                << endl;      // write the result from popular ones to minor ones
        if(haplotype_frequency<.01)
            maphl2rareFlag[mit->first]="r";
        else
            maphl2rareFlag[mit->first]="c";
    }
    for(i=0; i<nHAPs; ++i)
        haplotypeRC[i] = haplotypeRC[i] + " " + maphl2rareFlag[hapsl[i]];
    sst.clear();
    vector<string>().swap(hapsl);
    map<string, string>().swap(maphl2rareFlag);
    map<string, unsigned short>().swap(maphl2count);
    map<string, unsigned short>().swap(maphl2homo_count);
    stmp="";
}

void beagleOutput::clear()
{
//    animalIDs.clear();
    vector<string>().swap(animalIDs);
//    mapMKN2Line.clear();
    map<string, string>().swap(mapMKN2Line);
    vector<string>().swap(haplotypeRC);
}



/*
 *  generate results by given numberSNPs per window
 *
 *
 *
 */


void beagleOutput::generateResultS(geneticMap &gmp, unsigned numSNPperWin)
{
    static string outFile1, outFile2, outFile3, outFile4;
    string strchr;
    stringstream sst;
    sst << chr;
    sst >> strchr;
    sst.clear();

    // "hdout.<breed>.<chr>.totl", counts number of haplotype alleles in each window
    outFile1 = "hdout." + breed + "." + strchr + ".totl";
    ofstream ofstotl (outFile1.c_str());

    // "hdout.<breed>.<chr>.dtld", counts times of occurrance of each haplotype allele in each window
    outFile2 = "hdout." + breed + "." + strchr + ".dtld";
    ofstream ofsdtld (outFile2.c_str());

    // "hdout.<breed>.<chr>.racm", shows the structure of (ra)re / (c)o(m)mon haplotypes for each animal: 2 lines per animal.
    outFile3 = "hdout." + breed + "." + strchr + ".racm";
    ofstream ofsracm (outFile3.c_str());

    // "hdout.<breed>.<chr>.rprp", rare proportion
    outFile4 = "hdout." + breed + "." + strchr + ".rprp";
    ofstream ofsrprp (outFile4.c_str());

    if(!ofstotl || !ofsdtld || !ofsracm || !ofsrprp)
    {
        ofstotl.close(); ofsdtld.close(); ofsracm.close(); ofsrprp.close();
        cout << "XXXXXXXXXXXX cannot open files to write into." << endl;
        exit (-30);
    }
    unsigned            i, curWIN (0);
    string              theLine;
    //  go through the map (within the chr)
    map<unsigned long, string>::iterator mit;
    vector<string> queue;

    ofstotl << left
    << setw(23) << "number_of_diff_haps"
    << setw(16) << "prop(noh/allh)"
    << setw(20) << "mean_of_hap_occs"
    << setw(8)  << "chr"
    << setw(8)  << "window"
    << setw(8)  << "breed"
    << endl;
    ofsdtld << left
    << setw(19) << "hap_length"
    << setw(numSNPperWin+5) << "haplotype"
    << setw(8)  << "occr"
    << setw(16)  << "freq"
    << setw(8)  << "chr"
    << setw(8)  << "window"
    << setw(8)  << "breed"
    << setw(10)  << "homo_occr"
    << setw(10)  << "homo_freq"
    << endl;

    map<unsigned long, string> *theMap (&(gmp.mapPos2Name));
    i=0;
    for(mit = theMap->begin(); mit != theMap->end(); ++mit)
    {
//        curPOS  = mit->first;               // first  = pos in map
        theLine = mapMKN2Line[mit->second]; // second = markerName mapped with pos, theLine = the line of beagle output containing given markerName
        queue.push_back(theLine);
        if(++i==numSNPperWin)
        {
            if(queue.size()==0)
                ofstotl << left
                << setw(23) << 0
                << setw(16) << 0
                << setw(20) << 0
                << setw(8)  << chr
                << setw(8)  << curWIN
                << setw(8)  << breed
                << endl;
            else
                eligantCalcAndWrite(ofstotl, ofsdtld, queue, curWIN, numSNPperWin+5);
            vector<string>().swap(queue);
            // 20150603: fill the gap between two windows: haplotype starting from last SNP of previous window
            queue.push_back(theLine);
            i=0;
            curWIN++;
        }
    }
    if(queue.size()!=0)
    {
        eligantCalcAndWrite(ofstotl, ofsdtld, queue, curWIN, numSNPperWin+5);
        vector<string>().swap(queue);
    }
    cout << "                + output file [" << outFile1.c_str() << "] generated." << endl;
    cout << "                + output file [" << outFile2.c_str() << "] generated." << endl;

    unsigned r_count(0), c_count(0), j;
    for(unsigned i=0; i<haplotypeRC.size(); ++i) // write out the haplotype rare/common info
    {
        ofsracm << left << setw(24) << animalIDs[i/2] << haplotypeRC[i] << endl;
        for(j=0; j<haplotypeRC[i].length(); ++j)
        {
            if(haplotypeRC[i][j] == 'r')
                ++r_count;
            else if(haplotypeRC[i][j] == 'c')
                ++c_count;
        }
        if(i && i%2==0)
        {
            ofsrprp << left << setw(8) << breed << left << setw(8) << strchr << setw(8) << r_count/float(r_count+c_count) << endl;
            c_count = r_count = 0;
        }
    }
    cout << "                + output file [" << outFile3.c_str() << "] generated." << endl;
    cout << "                + output file [" << outFile4.c_str() << "] generated." << endl;
    //    cout << endl;
    ofstotl.close();
    ofsdtld.close();
    ofsracm.close();
    ofsrprp.close();
    //    animalIDs.clear();
    vector<string>().swap(animalIDs);
    //    mapMKN2Line.clear();
    map<string, string>().swap(mapMKN2Line);
    vector<string>().swap(haplotypeRC);
}
