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
//  beagleOutput.h
//  haplotype_diversity
//
//  Created by Hailin SU on 9/11/13.
//

#ifndef __haplotype_diversity__beagleOutput__
#define __haplotype_diversity__beagleOutput__

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include "geneticMap.h"
#include <fstream>

using namespace std;

class beagleOutput
{
public:
    short chr;                          // which chr
    vector<string>      animalIDs;
    map<string, string> mapMKN2Line;
    vector<string>      haplotypeRC;    // haplotype is rare/common

    string breed;

    void    load(string fileName, justAmap &jmp);
    string  getAnimalIDbyIndex(unsigned int);
    void    generateResult(geneticMap &gmp);
    void    generateResultS(geneticMap &gmp, unsigned numSNPperWin);
    void    generateResultW(geneticMap &gmp, float winSizeM);
    void    eligantCalcAndWrite(ofstream &ofstotl, ofstream &ofsdtld, vector<string> &queue, unsigned int curWIN, unsigned hapWidth);
    void    clear();
};

#endif /* defined(__haplotype_diversity__beagleOutput__) */
