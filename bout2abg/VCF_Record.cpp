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
//  VCF_Record.cpp
//  libHailins
//
//  Created by Hailin SU on 2/16/14.
//  Copyright (c) 2014 iastate. All rights reserved.
//

#include "VCF_Record.h"
#include <iomanip>

#define L(...)              lg.generalLog("I", __VA_ARGS__)

using namespace std;

namespace libcbk
{
    unsigned VCF_File::WINDOW_SIZE(0);
    const unsigned VCF_File::M1MB(1000000);  // 1,000,000

    void VCF_Record::clear()
    {
        chr=pos=0;
        ID="";
        vector<bool>().swap(vec_b_GenotypeLine);
    }

    void VCF_Record::loadFromSS(splitablestring ss)
    {
        vector<bool>().swap(vec_b_GenotypeLine);
        static vector<string> vecStrTmp;
        static stringstream    sst;
        static unsigned        i;

        vecStrTmp=ss.split('\t');
        // tab-splited
        //  0       1       2       3       4       5       6       7       8       9+
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BSHUSAM000004054458
        sst.str("");
        sst << vecStrTmp[0];
        sst >> chr;
        sst.clear();

        sst << vecStrTmp[1];
        sst >> pos;
        sst.clear();

        ID = vecStrTmp[2];
        for(i=9; i<vecStrTmp.size(); ++i)
        {
            // 0|1:1:0,1,0, [0]='0', [2]='1'
            vec_b_GenotypeLine.push_back(vecStrTmp[i][0]-'0');
            vec_b_GenotypeLine.push_back(vecStrTmp[i][2]-'0');
        }
        vector<string>().swap(vecStrTmp);
    }

    bool VCF_File::open()
    {
        bool ret(false);
        lastVCFR.chr=lastVCFR.pos=0;
        infs.open(fileName.c_str());
        if(infs)
        {
            FILE_SIZE_1p=0.01*tellFileSize(fileName);
            splitablestring aLine;
            unsigned linesRead(0);      // break after reading >200 lines
            while(getline(infs, aLine)) // skip header and read sampleIDs
            {
                if(aLine.size()<2)  // make sure I can read aLine[1]
                    continue;
                else if(aLine[0]=='#' && aLine[1]=='#')  // skip meta header
                    continue;
                else if(aLine[0]=='#')  // header, read and break
                {
                    // at least 10 columns should go
                    sx=aLine.split('\t');
                    if(sx.size()<10)
                    {
                        infs.close();
                        break;              // return false
                    }
                    for(unsigned i=9; i<sx.size(); ++i)
                        vecSampleIDs.push_back(sx[i]);
                    breed=vecSampleIDs[0].substr(0, 3);
                    ret=true;
                    break;
                }
                if(++linesRead>200)         // prevent large non-VCF file
                    break;
            }
        }
        return ret;
    }

    CQueue VCF_File::getNextQueue()
    {
        return getNextQueueF(1);   // by default, 1Mbp window
    }

    CQueue VCF_File::getNextQueueF(float windowSizeMbp)
    {
        if(WINDOW_SIZE==0)
            WINDOW_SIZE=windowSizeMbp*M1MB;// 1, 000, 000
        queue.clear();

        if(lastVCFR.chr==lastVCFR.pos && lastVCFR.chr==0)   // need to get 1st window
        {
            if(getline(infs, aLine))    // read 1st record, get current window
                lastVCFR.loadFromSS(aLine);
            else                        // not even could I read a record, just return empty
                return CQueue();
        }

        if(queue.size() && queue[0].chr!=lastVCFR.chr)  // different chr, drop the last SNP from previous window
            vector<VCF_Record>().swap(queue);
        queue.curChr=lastVCFR.chr;
        queue.curWindow=unsigned(lastVCFR.pos/WINDOW_SIZE);
        queue.push_back(lastVCFR);
        lastVCFR.clear();

        bool eof(true);
        while(getline(infs, aLine))
        {
            lastVCFR.loadFromSS(aLine);
            if( queue.curChr==lastVCFR.chr &&
               queue.curWindow==unsigned(lastVCFR.pos/WINDOW_SIZE))    // still same window
                queue.push_back(lastVCFR);
            else
            {
                eof=false;
                break;
            }
        }
        if(eof)
            lastVCFR.clear();
        curTellG=infs.tellg();
        return queue;
    }

    CQueue VCF_File::getNextQueueU(unsigned nSNPperWindow)
    {
        static unsigned short lastChr;
        static unsigned short lastWin;
        lastChr=queue.curChr;
        lastWin=queue.curWindow;
        queue.clear();
        static unsigned short recordsRead;
        if(getline(infs, aLine))    // read 1st record, get current chr
        {
            lastVCFR.loadFromSS(aLine);
            queue.curChr=lastVCFR.chr;
            if(queue.curChr==lastChr)
                queue.curWindow=lastWin+1;
            else
                queue.curWindow=0;
            queue.push_back(lastVCFR);
            recordsRead=queue.size();

            while(++recordsRead <= nSNPperWindow && getline(infs, aLine))
            {
                lastVCFR.loadFromSS(aLine);
                if(queue.curChr==lastVCFR.chr)    // still same chr
                    queue.push_back(lastVCFR);
                else
                    break;
            }
        }
        else
            return CQueue();    // could not read any more, return empty
        curTellG=infs.tellg();
        return queue;
    }

    void VCF_File::close()
    {
        queue.clear();
        vector<string>().swap(vecSampleIDs);
        vector<string>().swap(sx);
        infs.close();
    }
/*
    void VCF_File::hdhdat(string fn)
    {
        vb=true;
        fileName=fn;
        if(!open())
        {
            E("ss", "Unable to open VCF file to read:", fn.c_str());
            return;
        }
        L("sss", 1, "File", fileName.c_str(), "opened successfully.");
        L("ss", 2, "Breed         :", breed.c_str());
        L("si", 2, "NumberSamples :", vecSampleIDs.size());

        string outFile1, outFile2, outFile3, outFile4;

            // "hdout.<breed>.<chr>.totl", counts number of haplotype alleles in each window
        outFile1 = "hdout." + breed + ".totl";
        ofstream ofstotl (outFile1.c_str());

            // "hdout.<breed>.<chr>.dtld", counts times of occurrance of each haplotype allele in each window
        outFile2 = "hdout." + breed + ".dtld";
        ofstream ofsdtld (outFile2.c_str());

            // "hdout.<breed>.<chr>.racm", shows the structure of (ra)re / (c)o(m)mon haplotypes for each animal: 2 lines per animal.
        outFile3 = "hdout." + breed + ".racm";
        ofstream ofsracm (outFile3.c_str());

            // "hdout.<breed>.<chr>.rprp", rare proportion
        outFile4 = "hdout." + breed + ".rprp";
        ofstream ofsrprp (outFile4.c_str());

        if(!ofstotl || !ofsdtld || !ofsracm || !ofsrprp)
        {
            ofstotl.close(); ofsdtld.close(); ofsracm.close(); ofsrprp.close();
            cout << "XXXXXXXXXXXX cannot open files to write into." << endl;
            exit (-30);
        }

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
            << setw(38) << "haplotype"
            << setw(8)  << "occr"
            << setw(16) << "freq"
            << setw(8)  << "chr"
            << setw(8)  << "window"
            << setw(8)  << "breed"
            << setw(10) << "homo_occr"
            << setw(10) << "homo_freq"
            << endl;

        CQueue cq;
        unsigned i, j, homo_count;
        unsigned nSNPs(0);
        const unsigned long nHAPs(vecSampleIDs.size()*2);
        string stmp;
        vector<string> hapsl;    // string supports no limit
        vector<string> haplotypeRC;
        haplotypeRC.resize(nHAPs);     //denotes haplotype is rare/common
        for(i=0; i<nHAPs; ++i)
            haplotypeRC[i]="";
        map<string, unsigned short> maphl2count;
        map<string, unsigned short>::iterator mit;
        float meanx, haplotype_frequency, homo_freq;
        map<string, string> maphl2rareFlag;
        map<string, unsigned short> maphl2homo_count;
        unsigned lastp(0);
        unsigned char ctmp;
        do
        {
            cq=getNextQueue();
            if(cq.size()==0)
                break;
            if(lastp!=unsigned(curTellG/FILE_SIZE_1p))
            {
                lastp=unsigned(curTellG/FILE_SIZE_1p);
                progShow('P', lastp);
            }
            // *                      A1  A1  A2  A2  A3  A3 ...
             *  queue[0]    SNP1
             *  queue[1]    SNP2
             *  queue[2]    SNP3
             *  ...         ...
             * /
            //  form the haplotypes
            nSNPs = unsigned(queue.size());   // numberSNPs in the current window
            // make all haplotype strings the same length
            stmp.resize(nSNPs, ' ');
//            stmp[nSNPs]=0;
            for(i=0; i<nHAPs; ++i)
            {
                stmp="";
                ctmp=0;
                for(j=0; j<nSNPs; ++j)
                {
                    ctmp<<=1;
                    ctmp+=(unsigned char)(cq[j].vec_b_GenotypeLine[i]);
                    if((1+j)%8==0)
                    {
                        stmp.push_back(ctmp);
                        ctmp=0;
                    }
                }
                if(ctmp!=0)
                {
                    ctmp<<=(8-(j%8));
                    stmp.push_back(ctmp);
                }
                hapsl.push_back(stmp);
            }
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
                        else                                                        // existed
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
            << setw(8)  << cq.curChr
            << setw(8)  << cq.curWindow
            << setw(8)  << breed
            << endl;

            for (mit=maphl2count.begin(); mit!=maphl2count.end(); ++mit)
            {
                ofsdtld.setf(ios::uppercase);
                ofsdtld << left << setw(8) << setfill(' ') << nSNPs;
                ofsdtld.unsetf(ios::uppercase);

                stmp=mit->first;    // the hap string
                if(maphl2homo_count.find(stmp)==maphl2homo_count.end())                // not in the homo map
                    homo_count = 0;
                else
                    homo_count = maphl2homo_count[stmp];
                homo_freq = 2.0 * homo_count / nHAPs;

                haplotype_frequency = mit->second/(float)nHAPs;

                ofsdtld << setw(5+unsigned(nSNPs)) << decode(stmp, nSNPs)  // the hap string
                << setw(8)  << mit->second                  // count
                << setw(16) << haplotype_frequency
                << setw(8)  << cq.curChr
                << setw(8)  << cq.curWindow
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
            vector<string>().swap(hapsl);
            map<string, string>().swap(maphl2rareFlag);
            map<string, unsigned short>().swap(maphl2count);
            map<string, unsigned short>().swap(maphl2homo_count);
            stmp="";
        } while (true);
        progClear();
        cout << endl
        << "                + output file [" << outFile1.c_str() << "] generated." << endl
        << "                + output file [" << outFile2.c_str() << "] generated." << endl;

        unsigned r_count(0), c_count(0);
        ofsrprp << left << setw(24) << "AnimalID" << left << setw(8) << "rare/(rare+comm)" << endl;
        for(i=0; i<nHAPs; ++i) // write out the haplotype rare/common info
        {
            ofsracm << left << setw(24) << vecSampleIDs[i/2] << haplotypeRC[i] << endl;
            for(j=0; j<haplotypeRC[i].length(); ++j)
            {
                if(haplotypeRC[i][j] == 'r')
                    ++r_count;
                else if(haplotypeRC[i][j] == 'c')
                    ++c_count;
            }
            if(i%2==1)
            {
                ofsrprp << left << setw(24) << vecSampleIDs[i/2] << left << setw(8) << r_count/float(r_count+c_count) << endl;
                c_count = r_count = 0;
            }
        }
        cout << endl
        << "                + output file [" << outFile3.c_str() << "] generated." << endl
        << "                + output file [" << outFile4.c_str() << "] generated." << endl << endl;
        //    cout << endl;
        ofstotl.close();
        ofsdtld.close();
        ofsracm.close();
        ofsrprp.close();
        //    animalIDs.clear();

        close();
    }
// *
    string VCF_File::decode(string cs, unsigned len)
    {
        static string ret;
        static unsigned i;
        ret.resize(len, ' ');
        for(i=0; i<len; ++i)
            ret[i] = (cs[i/8] << (i%8)) & 128; // 0b 1000 0000
        return ret;
    }
*/
    void VCF_File::hdhdat(const string fn, const string type, const string val)
    {
        // vb=true;
        fileName=fn;
        if(!open())
        {
            E("ss", "Unable to open VCF file to read:", fn.c_str());
            return;
        }
        L("sss", 1, "File", fileName.c_str(), "opened successfully.");
        L("ss", 2, "Breed         :", breed.c_str());
        L("si", 2, "NumberSamples :", vecSampleIDs.size());

        string outFile1, outFile2, outFile3, outFile4;

        // "hdout.<breed>.<chr>.totl", counts number of haplotype alleles in each window
        outFile1 = "hdout." + breed + ".totl";
        ofstream ofstotl (outFile1.c_str());

        // "hdout.<breed>.<chr>.dtld", counts times of occurrance of each haplotype allele in each window
        outFile2 = "hdout." + breed + ".dtld";
        ofstream ofsdtld (outFile2.c_str());

        // "hdout.<breed>.<chr>.racm", shows the structure of (ra)re / (c)o(m)mon haplotypes for each animal: 2 lines per animal.
        outFile3 = "hdout." + breed + ".racm";
        ofstream ofsracm (outFile3.c_str());

        // "hdout.<breed>.<chr>.rprp", rare proportion
        outFile4 = "hdout." + breed + ".rprp";
        ofstream ofsrprp (outFile4.c_str());

        if(!ofstotl || !ofsdtld || !ofsracm || !ofsrprp)
        {
            ofstotl.close(); ofsdtld.close(); ofsracm.close(); ofsrprp.close();
            cout << "XXXXXXXXXXXX cannot open files to write into." << endl;
            exit (-30);
        }

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
        << setw(38) << "haplotype"
        << setw(8)  << "occr"
        << setw(16) << "freq"
        << setw(8)  << "chr"
        << setw(8)  << "window"
        << setw(8)  << "breed"
        << setw(10) << "homo_occr"
        << setw(10) << "homo_freq"
        << endl;

        CQueue cq;
        unsigned i, j, homo_count;
        unsigned long nSNPs(0);
        const unsigned long nHAPs(vecSampleIDs.size()*2);
        string stmp;
        vector<string> hapsl;    // string supports no limit
        vector<string> haplotypeRC;
        haplotypeRC.resize(nHAPs);     //denotes haplotype is rare/common
        for(i=0; i<nHAPs; ++i)
            haplotypeRC[i]="";
        map<string, unsigned short> maphl2count;
        map<string, unsigned short>::iterator mit;
        float meanx, haplotype_frequency, homo_freq;
        map<string, string> maphl2rareFlag;
        map<string, unsigned short> maphl2homo_count;
        unsigned lastp(0);
        stringstream sst;
        float par;
        sst << val;
        sst >> par;
        sst.clear();

        do
        {
            if(type=="-f")
                cq=getNextQueueF(par);
            else
                cq=getNextQueueU(unsigned(par));
            if(cq.size()==0)
                break;
            if(lastp!=unsigned(curTellG/FILE_SIZE_1p))
            {
                lastp=unsigned(curTellG/FILE_SIZE_1p);
                progShow('P', lastp);
            }
            /*                      A1  A1  A2  A2  A3  A3 ...
             *  queue[0]    SNP1
             *  queue[1]    SNP2
             *  queue[2]    SNP3
             *  ...         ...
             */
            //  form the haplotypes
            nSNPs = queue.size();   // numberSNPs in the current window
            // make all haplotype strings the same length
            stmp.resize(nSNPs, ' ');
            //            stmp[nSNPs]=0;
            for(i=0; i<nHAPs; ++i)
            {
                for(j=0; j<nSNPs; ++j)
                    stmp[j]=(unsigned char)(cq[j].vec_b_GenotypeLine[i]) + '0';
                hapsl.push_back(stmp);
            }
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
                        else                                                        // existed
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
            << setw(8)  << cq.curChr
            << setw(8)  << cq.curWindow
            << setw(8)  << breed
            << endl;

            for (mit=maphl2count.begin(); mit!=maphl2count.end(); ++mit)
            {
                ofsdtld.setf(ios::uppercase);
                ofsdtld << left << setw(8) << setfill(' ') << nSNPs;
                ofsdtld.unsetf(ios::uppercase);

                stmp=mit->first;    // the hap string
                if(maphl2homo_count.find(stmp)==maphl2homo_count.end())                // not in the homo map
                    homo_count = 0;
                else
                    homo_count = maphl2homo_count[stmp];
                homo_freq = 2.0 * homo_count / nHAPs;

                haplotype_frequency = mit->second/(float)nHAPs;

                ofsdtld << setw(5+unsigned(nSNPs)) << stmp  // the hap string
                << setw(8)  << mit->second                  // count
                << setw(16) << haplotype_frequency
                << setw(8)  << cq.curChr
                << setw(8)  << cq.curWindow
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
            vector<string>().swap(hapsl);
            map<string, string>().swap(maphl2rareFlag);
            map<string, unsigned short>().swap(maphl2count);
            map<string, unsigned short>().swap(maphl2homo_count);
        } while (true);
        progClear();
        cout << endl
        << "                + output file [" << outFile1.c_str() << "] generated." << endl
        << "                + output file [" << outFile2.c_str() << "] generated." << endl;

        unsigned r_count(0), c_count(0);
        ofsrprp << left << setw(24) << "AnimalID" << left << setw(8) << "rare/(rare+comm)" << endl;
        for(i=0; i<nHAPs; ++i) // write out the haplotype rare/common info
        {
            ofsracm << left << setw(24) << vecSampleIDs[i/2] << haplotypeRC[i] << endl;
            for(j=0; j<haplotypeRC[i].length(); ++j)
            {
                if(haplotypeRC[i][j] == 'r')
                    ++r_count;
                else if(haplotypeRC[i][j] == 'c')
                    ++c_count;
            }
            if(i%2==1)
            {
                ofsrprp << left << setw(24) << vecSampleIDs[i/2] << left << setw(8) << r_count/float(r_count+c_count) << endl;
                c_count = r_count = 0;
            }
        }
        cout << endl
        << "                + output file [" << outFile3.c_str() << "] generated." << endl
        << "                + output file [" << outFile4.c_str() << "] generated." << endl << endl;

        //    cout << endl;
        ofstotl.close();
        ofsdtld.close();
        ofsracm.close();
        ofsrprp.close();
        //    animalIDs.clear();
        
        close();
    }
}