import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg


   fs = [line.rstrip('\n') for line in open(inputfile)]
   f = [float(x) for x in fs]
   nx = len(fs)
   x = np.linspace(0,2*np.pi,nx)

   plt.plot(x,f,'b')
   plt.savefig(outputfile)
   plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])
