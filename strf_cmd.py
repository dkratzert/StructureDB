
import sys
import argparse

import time

from searcher import filecrawler, database_handler

#from searcher.spinner import Spinner


parser = argparse.ArgumentParser(description='Command line version of StructureFinder to collect cif files to a '
                                             'database.\n'
                                             'StructureFinder will search for cif files in the given directory(s) '
                                             'recursively.')
parser.add_argument("-d",
                    dest="dir",
                    metavar='"directory"',
                    type=str,
                    action='append',
                    help='Directory(s) where cif files are located.')
parser.add_argument("-e",
                    dest="ex",
                    metavar='"directory"',
                    type=str,
                    action='append',
                    help='Directory names to be excluded from the file search. Default is:\n'
                         '"ROOT", ".OLEX", "TMP", "TEMP", "Papierkorb", "Recycle.Bin" '
                         'Modifying -e option discards the default.')
parser.add_argument("-o",
                    dest="outfile",
                    metavar='"file name"',
                    type=str,
                    help='Name of the output database file. Default: "structuredb.sqlite"')


args = parser.parse_args()


ncifs = 0
try:
    if not args.dir:
        parser.print_help()
        sys.exit()
except IndexError:
    print("No valid search directory given.\n")
    print("Please run this as 'stdb_rmd [directory]'\n")
    print("stdb_cmd will search for .cif files in [directory] recoursively.")
else:
    dbfilename = 'structuredb.sqlite'
    if args.outfile:
        dbfilename = args.outfile
    # the command line version
    db = database_handler.DatabaseRequest(dbfilename)
    db.initialize_db()
    lastid = db.get_lastrowid()
    if not lastid:
        lastid = 0
    structures = database_handler.StructureTable(dbfilename)
    try:
        time1 = time.clock()
        for p in args.dir:
                ncifs = filecrawler.put_cifs_in_db(searchpath=p, excludes=args.ex,
                                                   structures=structures, lastid=lastid)
        time2 = time.clock()
        diff = time2 - time1
        m, s = divmod(diff, 60)
        h, m = divmod(m, 60)
        tmessage = '\nTotal {3} cif files. Duration: {0:>2d} h, {1:>2d} m, {2:>3.2f} s'
        print(tmessage.format(int(h), int(m), s, ncifs))
        print("Written to '{}'".format(dbfilename))
    except OSError as e:
        print("Unable to collect files:")
        print(e)
    #spinner.stop()
    pass