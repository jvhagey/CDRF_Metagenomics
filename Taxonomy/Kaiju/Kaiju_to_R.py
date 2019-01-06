#!/usr/bin/env python

## Jill Hagey
## University of California, Davis
## jvhagey@gmail.com
## https://github.com/jvhagey/
## 2018
## Given a Kaiju .summary file script will give a .csv file to read into R

#importing packages
import glob,os
import pandas as pd
import re
import time
import sys
import logging
import logging.handlers
from argparse import ArgumentParser

#Defining function to process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="exsum.py",description="Parses Kaiju .summary files into an R readable .csv file")
    parser.add_argument("-T",'--Taxonomic-level', dest="taxlevel", action="store", required=True, help='Pick One: phylum, class, order, family, genus, species (required)')
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None, required=True, help="Input directory name were .summary files are found (required)")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None, required=True, help="Output directory (required)")
    parser.add_argument("-f", "--force", dest="force", action="store_true", default=False, help="Force file overwriting")
    parser.add_argument("--nooverwrite", dest="nooverwrite", action="store_true", default=False, help="Don't override existing files")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None, help="Logfile location")
    return parser.parse_args()

#Defining function to create output directory if it doesn't exist
def make_outdir():
    """Make the output directory, if required.

    If the output directory already exists by default the script stops with an error.
    If you force the program to go on, you can either overwrite the existing directory, or not.
    The options turn out as the following, if the directory exists:

    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOOVERWRITE+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would overwrite existing " +
                         "files (exiting)", args.outdirname)
            sys.exit(1)
        elif args.noclobber:
            logger.warning("NOOVERWRITE: not actually deleting directory %s",
                           args.outdirname)
        else:
            logger.info("Removing directory %s and everything below it",
                        args.outdirname)
            shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        os.makedirs(args.outdirname)   #Make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOOVERWRITE+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)

if __name__ == '__main__':
    #Parse command-line
    args = parse_cmdline()
    #Set up logging
    logger = logging.getLogger('exsum.py: %s' %
                               time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    #Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging",
                         args.logfile)
            sys.exit(1)
    #Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s", args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)

    #getting name of file from taxlevel
    line = '*_kaiju_out_.summary'
    index = line.find('.summary')
    output_line = line[:index] + args.taxlevel + line[index:]
    #Calling in files
    os.chdir(args.indirname)
    infiles = glob.glob(output_line)
    #saving new header to be used later in rewriting file
    new_header = ("Percent_Reads\tReads\t") + args.taxlevel.capitalize()
    df_header = ("Sample Percent_Reads Reads") + args.taxlevel.capitalize()
    summary = [] #making blank dataframe to store counts of each gene from output of for loop
    for file in infiles:
        with open(file, "r+") as f, open(f.name + '.csv', 'w') as f2: #open file for reading and writing and antoher to write results into
            lines = f.readlines()
            f2.writelines(lines[2:-5]) #write lines from orignal file and write into new file removing the first 3 and last 5
        with open(f.name + '.csv', 'r') as f2:
            filedata = f2.read()
        with open(f.name + '.csv', 'w') as f3: #Read in the new file that was just written
            #Replacing spaces between two works and replace with _
            filedata = re.sub(r'([A-Za-z]) ([A-Za-z])', r'\1_\2',filedata)
            filedata = re.sub(" " , "", filedata)
            f3.write(new_header) #writing cleaned lines into file
            f3.write('\n')
            f3.writelines(filedata)
        filename = f.name + '.csv'
        df = pd.read_csv(filename, sep='\t', header=0, usecols=["Percent_Reads", "Reads", args.taxlevel.capitalize()])
        df['Sample'] = f.name
        #editing output line name to remove *
        output_line2 = output_line.replace('*', '')
        df['Sample'] = df['Sample'].str.replace('_S(.*?)' + output_line2, '') #editing "Sample_Names"
        df['Taxanomic_Level'] = f.name
        df['Taxanomic_Level'] = args.taxlevel.capitalize()
        #Numbers at beginning of file are Month-Farm-Sample
        df['Farm'] = pd.np.where(df.Sample.str.contains("8"), "Farm_8", #Adding column that will have the farm name
                           pd.np.where(df.Sample.str.contains("6"), "Farm_6",
                           pd.np.where(df.Sample.str.contains("5"), "Farm_5",
                           pd.np.where(df.Sample.str.contains("1"), "Farm_1", "No_name"))))
        summary.append(df) #store DataFrame in list
    summary = pd.concat(summary)
#write DataFrame to comma separated file (.csv) with file name and TIGRFAM counts
os.chdir(args.outdirname)
name = args.taxlevel + '_summary.csv'
summary.to_csv(name, sep=',' ,index=False)

#Report that we've finished
logger.info("Done: %s.", time.asctime())
logger.info("Time taken: %.2fs", (time.time() - t0))
