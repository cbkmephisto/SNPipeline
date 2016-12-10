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
//  logger.h
//  simple logger, using just std::out
//
//  Created by Hailin Su on Wed, Aug 28, 2013.
//

/*	comment out these two lines to stop showing the I() log */
#ifndef LOGGER_H
#define LOGGER_H

#include<fstream>
#include<sstream>
using namespace std;

namespace libcbk
{
    class Logger
    {
    private:
        clock_t         startTime;
        time_t          t;
        tm              *lt;
        va_list         argp;	/*	structure storing param list	*/
        char*           para;	/*	temp pointer for current param	*/
        unsigned short  i;
        char*           szTypes;
        union           Printable_t
        {
            long long   i;
            double      f;
            char*       s;
        } Printable;
        bool            outStrI, outStrE, outStrF;
        char buff[64];

    public:
        ofstream filex;

        void generalLog(string outStr, string typeStr, int indLvl, ... );
        void initLogger();
        void finalizeLogger();
        void initF(string logFileName);
        void closF();
    };
}
#endif
