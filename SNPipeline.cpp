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
//  SNPipeline.cpp
//  libHailins
//
//  Created by Hailin SU on 10/3/13.
//
//  SNPipeline - step 1 - creating ab-genotype files

#include "appInfo.h"
#include "data.h"
//#include "finalReport.h"
//#include "affyInput.h"


using namespace std;
namespace libcbk
{
    class SNPipeline : public appInfo
    {
    public:
        SNPipeline(const int argc, const char * argv[]);
    private:
        void theMain(const int argc, const char * argv[]);
    };

    SNPipeline::SNPipeline(const int argc, const char * argv[])
    {
        progName    = "SNPipeline";
        version     = "2015-08-13 11:27 CDT";
        created     = "2013-10-01";
        updated     = "2015-08-12";
        descrip     = "SNPipeline - step 1. Creates AB-genotype file from finalReport/affyInput file(s). An optional space/tab-delimited >=2-col file named 'xref' is considered: [0]->[last].\n\n";
        descrip    += "[Warnings]\n";
        descrip    += "    - ALWAYS MAKE THE INPUT FILES INTO DIFFERENT FOLDERS!! If 2 or more finalReport/affyInput files with the same [CONTENT] and numberSNPs exist in the SAME DIRectory, the resulting ab-genotype file for these files would have the same fileName, and only the last processed file will be in the folder-related block of the resulting ab-genotype file. It's BLOCK-REPLACING, neither APPENDING nor MERGING!\n";
        descrip    += "\n[xref]\n";
        descrip    += "    - A space/tab-delimited->=2column-per-line sampleID cross-reference file in the same directory with finalReport/affyInput file named 'xref' will be used if exists.\n";
        descrip    += "    - This program will replace all the spaces in the sampleIDs in the input finalReport/affyInput files by underscore char('_'), thus no white space char(' ') is expected in either column of the IDs in the xref file. Actually ' ' is a delimit char, while tab('\t', '\\t') is the other.\n";
        descrip    += "\n[Progress % Legend]\n";
        descrip    += "    C for  pre-checking\n";
        descrip    += "    P for post-checking\n";
        descrip    += "    L for loading genotype\n";
        descrip    += "    T for creating temp file\n";
        descrip    += "    W for writing\n";

//        displayInfo();
//        lg.initLogger();
        theMain(argc, argv);
//        lg.finalizeLogger();
        /*
         updates
         2013-11-13 make seperated class: finalReport.cpp
         2013-11-20 space-delimitered 2-col xref file considered
         2013-12-17 - duplicated SampleIDs allowed after cross referencing
                    - reprocess any finalReport file if processed finalReport files were given as parameter
                    - do not put more than 1 finalReport file in the same folder: previously processed one would be discarded.
         2014-01-23 - added pipeline for affymetrix genotype format
                    - xref can be space/tab delimited
                    - a warning message will be printed on the screen if generating a NEW ab-genotype file
                    - a character denoting what the program is doing was added to the progress %
                    - timing for jobs
         2014-01-29 - class 'finalReport' and 'affyInput' were written into one class 'data' due to many duplicate functions
                    - finalReport speed up ~ 40% (135s -> 78s for instance)
         2014-02-06 fixed a potentially index out of boundary problem for not checking the nSNPs and nSamples in the header part
         2015-08-13 - duplicated sampleIDs are allowed
                    - fixed problems caused by too-long-sampleIDs
                    - xref could be more than 2 columns, xrefing [0]->[last]
        */
    }

    void SNPipeline::theMain(const int argc, const char * argv[])
    {
        makeSure(argc >= 2, argv[0], "'Red Angus 50K 29feb2012_FinalReport.txt' (*.txt)");
        for(unsigned argind=1; argind<argc; ++argind)
        {
//            I("s", 1, "******************************************************");
            I("sss", 1, "{", argv[argind], "}");
            static Data dt;
            dt.pipeline(argv[argind]);
            /*
            int ret;

            if(isIllumina(argv[argind]))                                        // Illumina chips
            {
                static finalReport fr;
                fr.clear();

                ret=fr.preCheck(argv[argind]);
                //  > 0, go ahead
                //  < 0, ERROR
                if(ret<0)
                    break;

                ret=fr.postCheck();
                //  == 0, info already exists, continue
                //   < 0, ERROR
                if(ret<0)
                    break;
                else if(ret==1) // simple append
                    fr.FR2AB_simpleAppend();
                else            // merging append
                    fr.FR2AB_mergingAppend();
                fr.clear();
            }
            else                                                                // Affymetrix
            {
                static affyInput af;
                af.clear();

                ret=af.preCheck(argv[argind]);
                //  > 0, go ahead
                //  < 0, ERROR
                if(ret<0)
                    break;

                ret=af.postCheck();
                //  == 0, info already exists, continue
                //   < 0, ERROR
                if(ret<0)
                    break;
                else
                {
//                    af.load();  // load data into memory
                    if(ret==1)  // simple append
                        af.AF2AB_simpleAppend();
                    else        // merging append
                        af.AF2AB_mergingAppend();
                }
                af.clear();
            }
            */
        }
    }
}

int main(int argc, const char * argv[])
{
    libcbk::SNPipeline(argc, argv);
    return 0;
}
