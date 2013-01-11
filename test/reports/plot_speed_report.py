

import sys
import numpy

from optparse import OptionParser

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def select(row, cols_spec):
    subspecs = cols_spec.split(',')
    subspecs = map(str.strip, subspecs)
    cols = []
    for subspec in subspecs:
        parts = subspec.split('-')
        if len(parts) == 1:
            cols.append(int(parts[0]))
        else:
            if len(parts[1]) == 0:
                end = len(row)
            else:
                end = int(parts[1])
                if end < 0:
                    end = len(row) + end
            cols.extend(range(int(parts[0]), end))
    for index in cols:
        yield row[index]
        
def plot_speed_report(input_filename = None, output_filename = None, cols = '0-'):
    with open(input_filename, 'r') as stream:
        lines = stream.readlines()
    header = None
    x = []
    data = []
    for line in lines:
        if line.startswith('#'):
            header_for_next_line = line[1:].split(',')
            header_for_next_line = list(select(header_for_next_line[2:], cols))
            if not header is None:
                if not header == header_for_next_line:
                    raise Exception("data does not have same header")
            header = header_for_next_line
        else:
            parts = line.split(',')
            x.append(int(parts[0]))
            numbers = map(lambda x : float(x), parts[2:])
            data.append(list(select(numbers, cols)))
    x = numpy.asarray(x)
    data = numpy.asarray(data)
    print data.shape
    
    figure = pyplot.figure(figsize=(9, 4))
    subplot = pyplot.subplot(1,2,1)
    
    handles = subplot.plot(x,data)

    subplot.legend(
        handles, 
        header,
        loc='center left', 
        bbox_to_anchor=(1.05, 0.5),
        ncol=1, 
        fancybox=False, 
        shadow=False)
    pyplot.loglog()
    if output_filename is None:
        pyplot.show()
    else:
        pyplot.savefig(output_filename)
                       
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-o", 
        default = None,
        dest="output_filename",
        help="save figure to output, by default it will display it",
        type="string"
    )
    result.add_option(
        "-i", 
        default = 'report.csv',
        dest="input_filename",
        help="name of the file to load the data from",
        type="string"
    )
    result.add_option(
        "--cols", 
        default = '0-',
        dest="cols",
        help="columns to plot, can by 1,2,3 or 0-3 or 0-5, 6, 3",
        type="string"
    )
    
    return result
    
    
    
if __name__ == '__main__':
    options, arguments = new_option_parser().parse_args()
    plot_speed_report(**options.__dict__)
