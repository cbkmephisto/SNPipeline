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
//  logger.cpp
//  simple logger, using just std::out
//
//  Created by Hailin Su on Wed, Aug 28, 2013.
//

#include "logger.h"
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <cstdarg>
#include <ctime>

using namespace std;

namespace libcbk
{
    //  generalLog() takes a format string 'typeStr' of the form
    //   "ifcs", where each character specifies the
    //   type of the argument in that position.
    //
    //  i = int
    //  f = float
    //  c = char    // no use
    //  s = string (char *)
    //
    //  Following the format specification is a variable
    //  list of arguments. Each argument corresponds to
    //  a format character in the format string to which
    //  the szTypes parameter points
    //
    //  A second parameter 'outStr' stands for the output stream
    //
    //      I = cout
    //      E = cerr
    //      F = filex
    //
    // EXAMPLE
    //
    // generalLog("issif", "I", 55, "String1", "String2", 66, 7.6f);
    //
    //      This statement will print "55 String1 String2 66 7.6"
    //          to stdout(I), stderr(E), or a log file(F). In this case it was 'I'
    //
    // ref: http://msdn.microsoft.com/en-us/library/fxhdxye9.aspx
    //      http://stackoverflow.com/questions/3530771/passing-variable-arguments-to-another-function-that-accepts-a-variable-argument
    /*
     use #define section listed below to make it more simple to use

     #define I(...)              lg->generalLog("I", __VA_ARGS__)
     #define E(str, ...)         lg->generalLog("E", str, 1, __VA_ARGS__)
     #define F(str, ...)         lg->generalLog("F", str, 1, __VA_ARGS__)
     */
    void Logger::generalLog(string outStr, string typeStr, int indLvl, ... )
    {
        outStrI=outStrE=outStrF=false;
        szTypes = (char*)outStr.c_str();
        for(i=0; szTypes[i]!='\0'; i++)
        {
            switch(szTypes[i])
            {
                case 'I': case 'i':
                    outStrI=true;
                    break;
                case 'E': case 'e':
                    outStrE=true;
                    break;
                case 'F': case 'f':
                    outStrF=true;
                default: ;
            }
        }

        t = time(NULL);
        lt = localtime(&t);

        if(outStrI)
            cout    << " + "
            << setfill('0') << setw(2) << lt->tm_hour << ":"
            << setfill('0') << setw(2) << lt->tm_min  << ":"
            << setfill('0') << setw(2) << lt->tm_sec << " | ";
        else if(outStrE)
            cerr    << " X "
            << setfill('0') << setw(2) << lt->tm_hour << ":"
            << setfill('0') << setw(2) << lt->tm_min  << ":"
            << setfill('0') << setw(2) << lt->tm_sec << " | ";
        else if(outStrF && filex)
            filex    << " + "
            << setfill('0') << setw(2) << lt->tm_hour << ":"
            << setfill('0') << setw(2) << lt->tm_min  << ":"
            << setfill('0') << setw(2) << lt->tm_sec << " | ";

        //	macro va_start, let argp point to the first param.
        //		noMsg is not a parameter in the list
        va_start(argp, indLvl);

        for(i=1; i<indLvl; i++)
            if(outStrI)
                cout    << "  ";
            else if(outStrE)
                cerr    << "  ";
            else if(outStrF && filex)
                filex   << "  ";

        if(indLvl>1)
        {
            if(outStrI)
                cout    << "- ";
            else if(outStrE)
                cerr    << "- ";
            else if(outStrF && filex)
                filex   << "- ";
        }
        szTypes = (char*)typeStr.c_str();
        // Step through the list.
        for(i = 0; szTypes[i] != '\0'; i++)
        {
            switch( szTypes[i] )
            {   // Type to expect.
                case 'i': case 'I':
                    Printable.i = va_arg( argp, int );
                    if(outStrI)
                        cout << Printable.i << " ";
                    else if(outStrE)
                        cerr << Printable.i << " ";
                    else if(outStrF && filex)
                        filex << Printable.i << " ";
                    break;

                case 'f': case 'F':
                    Printable.f = va_arg( argp, double );
                    if(outStrI)
                        cout << Printable.f << " ";
                    else if(outStrE)
                        cerr << Printable.f << " ";
                    else if(outStrF && filex)
                        filex << Printable.f << " ";
                    break;
                case 's': case 'S':
                    Printable.s = va_arg(argp, char*);
                    if(outStrI)
                        cout << Printable.s << " ";
                    else if(outStrE)
                        cerr << Printable.s << " ";
                    else if(outStrF && filex)
                        filex << Printable.s << " ";
                default:
                    ;
            }
        }
        va_end(argp);	/* set argp to NULL	*/
        if(outStrI)
            cout << endl;
        else if(outStrE)
            cerr << endl;
        else if(outStrF && filex)
            filex << endl;
    }
    
    void Logger::initLogger()
    {
        t = time(NULL);
        lt = localtime(&t);
        sprintf(buff, "Starting %d-%02d-%02d, with PID = %d",
                (1900+lt->tm_year),
                (lt->tm_mon+1),
                (lt->tm_mday),
                getpid());
        //    cout << buff;
        generalLog("I", "s", 1, buff);
        startTime = clock();
    }
    
    void Logger::finalizeLogger()
    {
        t = time(NULL);
        lt = localtime(&t);
        
        sprintf(buff, "Ending   %d-%02d-%02d, %lu seconds.", (1900+lt->tm_year), (lt->tm_mon+1), (lt->tm_mday), (clock()-startTime)/CLOCKS_PER_SEC);
        generalLog("I", "s", 1, buff);
    }

    void Logger::initF(string logFileName)
    {
        if(filex)
            filex.close();
        filex.open(logFileName.c_str());
    }

    void Logger::closF()
    {
        if(filex)
            filex.close();
    }
}