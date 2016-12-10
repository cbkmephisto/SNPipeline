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
//  SNP.h
//  SNPipeline
//
//  Created by Hailin SU on 3/24/14.
//

#ifndef __SNPipeline__SNP__
#define __SNPipeline__SNP__

namespace libcbk
{
    class SNP
    {
    public:
        string ID;
        unsigned chr;
        unsigned long pos;
        static bool SNPcomp(const SNP SNP1, const SNP SNP2);
        SNP(const string xID, const unsigned xchr, const unsigned long xpos)
        {
            this->ID=xID;
            this->chr=xchr;
            this->pos=xpos;
        }

        bool gt(const SNP SNPx)
        {
            if( chr>SNPx.chr || (pos>SNPx.pos && chr==SNPx.chr && ID != SNPx.ID))
                return true;
            return false;
        }

        bool identicalTo(const SNP SNPx)
        {
            if((pos==SNPx.pos && chr==SNPx.chr) || ID == SNPx.ID)
                return true;
            return false;
        }

        bool equalsbyID(const SNP SNPx)
        {
            if(ID == SNPx.ID)
                return true;
            return false;
        }

        bool equalsbyPos(const SNP SNPx)
        {
            if(pos==SNPx.pos && chr==SNPx.chr)
                return true;
            return false;
        }
    };

    bool SNP::SNPcomp(const SNP SNP1, const SNP SNP2)
    {
        //  same chromosomes compare map positions
        if (SNP1.chr == SNP2.chr)
            return SNP1.pos < SNP2.pos ? true: false;
        else
            return SNP1.chr < SNP2.chr ? true: false;
    }

}
#endif /* defined(__SNPipeline__SNP__) */
