#!/home/msiebert/anaconda2/bin/python2.7

import sys
import ntt
from optparse import OptionParser
import util
import quick_reduc


description = "> Fast reduction of efosc spectra "
usage = "%prog  \t raw_spectrum [option] "

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog " + str(ntt.__version__))
    parser.add_option("-d", "--dispersion", dest="dispersion", action="store_true",
                      help='chose interactively dispersion line')
    parser.add_option("-t", "--trace", dest="trace", action="store_true", help=' trace extraction with another frame')
    parser.add_option("-s", "--sens", dest="sens", default='', type="str", help=' use sensitivity curve from this list')
    parser.add_option("-a", "--arc", dest="arc", default='', type="str", help=' use arc from this list ')
    parser.add_option("-i", "--interactive", dest="interactive", action="store_true")
    parser.add_option("-c", "--cosmic", dest="cosmic", action="store_true")

    option, args = parser.parse_args()
    if len(args) < 1:
        sys.argv.append('--help')
    option, args = parser.parse_args()

    _cosmic = option.cosmic
    if not _cosmic:   _cosmic = False
    _interactive = option.interactive
    if option.trace == None:
        _trace = 'No'
    else:
        _trace = 'Yes'
    _dispersionline = option.dispersion
    if not _dispersionline:   _dispersionline = False
    if option.arc:
        _listarc = util.readlist(option.arc)
    else:
        _listarc = ''
    if option.sens:
        _listsens = util.readlist(option.sens)
    else:
        _listsens = ''
    # imglist = util.readlist(args)
    imglist = args
    outputfile = quick_reduc.reduce(imglist, _listsens, _listarc, _trace, _dispersionline, _cosmic,
                                                    _interactive)
    print '\n#######################################\n### end of reduction \n### output files:'
    for img in outputfile:
        print img


