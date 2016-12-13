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
//  fout2haplotype.h
//  libHailins
//
//  Created by Hailin SU on 4/21/14.
//

#ifndef __libHailins__fout2haplotype__
#define __libHailins__fout2haplotype__

#include "../appInfo.h"
#include <vector>

using namespace std;

namespace libcbk
{
    // ######################################################################### class HapWindow
    class HapWindow
    {
    public:
        unsigned indexStart;
        unsigned indexEnd;
        unsigned i;
        char     ctmp;

        string chr;
        string win;

//        string dipl;    // diplotype
        string hapP;    // 1st line
        string hapM;    // 2nd line

        HapWindow(unsigned st, unsigned en, const string c, const string w)
        : indexStart(st), indexEnd(en), chr(c), win(w) {}

        void load(string &str)
        {
//            dipl=(*str).substr(indexStart, indexEnd-indexStart+1);
            hapP=hapM=chr;
            hapP+="_";
            hapM+="_";

            hapP+=win;
            hapM+=win;

            hapP+="_";
            hapM+="_";

//            hapP=hapM="";
            for(i=indexStart; i<=indexEnd; ++i)
            {
                ctmp=str.at(i);
                /* valid: 0, 2, 3, 4, 1
                A1->0, A2->1, unphased: +-
                 0: A1A1    00
                 1: Unphased heterozygous
                 2: A2A2    11
                 3: A1A2    01
                 4: A2A1    10
                 5: missing
                 6: A1–
                 7: A2–
                 8: –A1
                 9: –A2
                 */
                switch (ctmp)
                {
                    case '0':
                        hapP+="0";
                        hapM+="0";
                        break;
                    case '1':
                        hapP+="+"; // +
                        hapM+="-"; // -
                        break;
                    case '2':
                        hapP+="1";
                        hapM+="1";
                        break;
                    case '3':
                        hapP+="0";
                        hapM+="1";
                        break;
                    case '4':
                        hapP+="1";
                        hapM+="0";
                        break;
                    case '5':
                        hapP+="X"; // missing
                        hapM+="X"; // missing
                        break;
                    default:
                        cerr << " *ERROR*: Invalid string parsed in HapWindow::load(), aborting." << endl;
                        cerr << "  ----->  " << str.substr(indexStart, indexEnd-indexStart+1) << endl;
                        exit(-10);
                }
            }
        }
    };

    // ######################################################################### class VecWindows
    class VecWindows: public vector<HapWindow>
    {
    public:
        bool            valid;
        unsigned        i;
        string          lineP;
        string          lineM;
        string          lineI;
        stringstream    sst;

        VecWindows(string snp_info_file_name)    // doing the init format
        {
            
            //            const unsigned W1M (1000000);   // 1M
            unsigned W1M (0);   // 1M
            ifstream infs ("windowSize");
            if(infs)
            {
                infs >> lineI;
                sst << lineI;
                sst >> W1M;
                lineI="";
                sst.clear();
                infs.close();
            }
            
            if(!W1M)
                W1M=1000000;    // default: 1M

            cout << "        using window size in bp: " << W1M << endl;
            valid=false;
            infs.open(snp_info_file_name.c_str());
            if(!infs)
            {
                cerr << " *ERROR* in VecWindows::init(): Invalid snp_info.txt file, aborting." << endl;
                cerr << "  ----->  " << snp_info_file_name << endl;
                exit (-20);
            }
            getline(infs, lineP);   // ignoring header line

            unsigned curPos, cb(W1M);
            unsigned lastChr(9999);
            unsigned curChr;
            unsigned curWin(0);

            unsigned mapIndex(0);
            unsigned mapStart(0);

            appInfo ai;
            unsigned long f1p (0.01*ai.tellFileSize(snp_info_file_name));
            if(f1p<1)
                f1p=1;
            ai.progShow('M', 0);
            while(infs >> lineP >> curChr >> curPos)//   ID  Chr pos
            {
                // get rid of the rest of the line
                getline(infs, lineM);

                // first read
                if(lastChr!=curChr && lastChr==9999)
                {
                    // update the new win
                    mapStart=mapIndex;
                    curWin=0;
                    cb=W1M;
                    lastChr=curChr;
                }
                // new chr
                else if(lastChr!=curChr && lastChr!=9999)
                {
                    // insert the last win
                    sst << lastChr;
                    sst >> lineP;   // chr
                    sst.clear();
                    sst << curWin;
                    sst >> lineM;   // win
                    sst.clear();
                    push_back(HapWindow(mapStart, mapIndex-1, lineP, lineM));
                    // update the new win
                    mapStart=mapIndex;
                    curWin=0;
                    cb=W1M;
                    lastChr=curChr;
                }
                // new win
                else if(cb<=curPos)
                {
                    // insert the last win
                    sst << curChr;
                    sst >> lineP;   // chr
                    sst.clear();
                    sst << curWin;
                    sst >> lineM;   // win
                    sst.clear();
                    push_back(HapWindow(mapStart, mapIndex-1, lineP, lineM));
                    // update the new win
                    // 20150604: mapStart-1: filling the gap between two adjacent windows
                    mapStart=mapIndex-1;
                    ++curWin;
                    cb=W1M*(curWin+1);
                }
                if(++mapIndex % 7000 == 1)
                    ai.progShow('M', infs.tellg()/f1p);
            }
            // push back the last window
            sst << curChr;
            sst >> lineP;   // chr
            sst.clear();
            sst << curWin;
            sst >> lineM;   // win
            sst.clear();
            push_back(HapWindow(mapStart, mapIndex-1, lineP, lineM));

            infs.close();
            ai.progClear();
            valid=true;
        }

        void load(string &pStr) // accept the whole line of genotype.txt
        {
            if(!valid)
            {
                cerr << " *ERROR* in VecWindows::load(): Invalid initialization before loading, aborting." << endl;
                exit (-21);
            }

            sst << pStr;
            sst >> lineM;           // store the animalID
            sst >> lineP >> lineI;  // >> chip >> genotype
            sst.str("");
            sst.clear();
            lineP=lineM;    // get the animalID

            for(i=0; i<size(); ++i)
            {
                (*this)[i].load(lineI);

                lineP+=" ";
                lineM+=" ";

                lineP+=(*this)[i].hapP;
                lineM+=(*this)[i].hapM;
            }
        }

        void load(string &pStr, string &c, string &w) // accept the whole line of genotype.txt
        {
            if(!valid)
            {
                cerr << " *ERROR* in VecWindows::load(): Invalid initialization before loading, aborting." << endl;
                exit (-21);
            }

            sst << pStr;
            sst >> lineM;           // store the animalID
            sst >> lineP >> lineI;  // >> chip >> genotype
            sst.str("");
            sst.clear();
            lineP=lineM;    // get the animalID

            for(i=0; i<size(); ++i)
            {
                if((*this)[i].win==w && (*this)[i].chr==c)
                {
                    (*this)[i].load(lineI);
                    
                    lineP+=" ";
                    lineM+=" ";

                    lineP+=(*this)[i].hapP;
                    lineM+=(*this)[i].hapM;
                    break;
                }
            }
        }
    };
}

#endif /* defined(__libHailins__fout2haplotype__) */
