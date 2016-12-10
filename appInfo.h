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
//  appInfo.h
//  libHailins
//
//  Created by Hailin SU on 9/20/13.
//

#ifndef __libHailins__appInfo__
#define __libHailins__appInfo__

#include <iostream>
#include <string>
#include "logger.h"
#include <cstdlib>

using namespace std;

namespace libcbk
{
    class appInfo
    {
    public:
        Logger      lg;
        string      progName;//    = "libHailins";
        string      version;//     = "init";
        string      created;//     = "N/A";
        string      updated;//     = "N/A";
        string      descrip;//     = "do something";

        char curp;
        const string pb, bl;
        unsigned xxi;

        appInfo() : pb ("....1....2....3....4....5....6....7....8....9....+"),
                    bl ("                                                                      "),
                    curp(-1){}

        void displayInfo(){}

        void progShow(const unsigned curPos)
        {
            if(curPos>=100) return;
            if(curp<0)
            {
                cout << "      x 10% :   0" << flush;
                curp=0;
            }
            for(xxi=curp; xxi<curPos/2; ++xxi)
                cout << pb[xxi] << flush;
            curp=xxi;
        }

        void progShow(const char ch, const unsigned curPos)
        {
            if(curPos>=100) return;
            if(curp<0)
            {
                cout << "    " << ch << " x 10% :   0" << flush;
                curp=0;
            }
            for(xxi=curp; xxi<curPos/2; ++xxi)
                cout << pb[xxi] << flush;
            curp=xxi;
        }

        void progClear()
        {
            if(curp>=0)
            {
                for(xxi=curp; xxi<50; ++xxi)
                    cout << pb[xxi];
            }
            cout << '\r' << bl << flush << '\r';
            curp=-1;
        }
/*
        void progShow(unsigned curPos)
        {
            cout
            << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
            << "        " << curPos << "%" << flush;
        }

        void progShow(char ch, unsigned curPos)
        {
            cout
            << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
            << "    " << ch << "    " << curPos << "%" << flush;
        }

        void progClear()
        {
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << flush;
        }
*/
        unsigned long long tellFileSize(string fileName)
        {
            ifstream xinfs (fileName.c_str());
            xinfs.seekg(0,std::ios::end);    // Seek to the end
            unsigned long long ret (xinfs.tellg());
            xinfs.close();
            return ret;
        }

        void usage(const string& argv0, const string& para)
        {
            cout << endl
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
            << "        " << argv0 << " " << para << endl << endl;
        }

        void makeSure(const bool testMe, const string& progN, const string& para)
        {
            if(!testMe)
            {
                usage(progN, para);
                exit(-101);
            }
        }
    };
}


#define I(...)              lg.generalLog("I", __VA_ARGS__)
#define E(str, ...)         lg.generalLog("E", str, 1, __VA_ARGS__)
#define F(str, ...)         lg.generalLog("F", str, 1, __VA_ARGS__)

#endif /* defined(__libHailins__appInfo__) */
