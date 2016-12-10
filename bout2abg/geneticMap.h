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
//  geneticMap.h
//  haplotype_diversity
//
//  Created by Hailin SU on 9/11/13.
//  Copyright (c) 2013 iastate. All rights reserved.
//

#ifndef __haplotype_diversity__geneticMap__
#define __haplotype_diversity__geneticMap__

#include <iostream>
#include <string>
#include <map>
#include <set>

using namespace std;

class justAmap
{
public:
    std::map<std::string, short> mapMKN2chr;                // markerName -> chr
    std::map<std::string, std::string> mapMKN2line;                // markerName -> line of map file
};

class geneticMap
{
public:
    unsigned short              chr;
    map<unsigned long, string>  mapPos2Name;  //  pos -> markerName

    static map<short, geneticMap> constructGeneticMap (string incf, string mapf, justAmap &jmp); //  chr -> geneticMap
};

#endif /* defined(__haplotype_diversity__geneticMap__) */
