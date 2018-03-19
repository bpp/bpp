#!/usr/bin/env python

from subprocess import Popen, PIPE, call

import sys, stat, os
import time

# define path to BPP binary

opt_bpp_bin = "$HOME/GIT/bpp/src/bpp"

# define each testbed and 

opt_testsuite_small_desc = "small"
opt_testsuite_small = [                 # [path-to-test,description]
   ["testbed/small/1", "small-A00-1"],
   ["testbed/small/2", "small-A00-2"],
   ["testbed/small/3", "small-A00-3"],
   ["testbed/small/4", "small-A00-4"],
   ["testbed/small/5", "small-A00-5"],
   ["testbed/small/6", "small-A00-6"],
   ["testbed/small/7", "small-A00-7"],
   ["testbed/small/8", "small-A00-8"],
   ["testbed/small/9", "small-A00-9"],
   ["testbed/small/10","small-A00-10"],
   ["testbed/small/11","small-A00-11"],
   ["testbed/small/12","small-A00-12"],
   ["testbed/small/13","small-A00-13"],
   ["testbed/small/14","small-A00-14"],
   ["testbed/small/15","small-A00-15"],
   ["testbed/small/16","small-A00-16"],
   ["testbed/small/17","small-A01-17"],
   ["testbed/small/18","small-A01-18"],
   ["testbed/small/19","small-A01-19"],
   ["testbed/small/20","small-A01-20"],
   ["testbed/small/21","small-A01-21"],
   ["testbed/small/22","small-A01-22"],
   ["testbed/small/23","small-A01-23"],
   ["testbed/small/24","small-A01-24"],
   ["testbed/small/25","small-A01-25"],
   ["testbed/small/26","small-A01-26"],
   ["testbed/small/27","small-A01-27"],
   ["testbed/small/28","small-A01-28"],
   ["testbed/small/29","small-A01-29"],
   ["testbed/small/30","small-A01-30"],
   ["testbed/small/31","small-A01-31"],
   ["testbed/small/32","small-A01-32"],
   ["testbed/small/33","small-A01-33"],
   ["testbed/small/34","small-A01-34"],
   ["testbed/small/35","small-A01-35"],
   ["testbed/small/36","small-A01-36"],
   ["testbed/small/37","small-A01-37"],
   ["testbed/small/38","small-A01-38"],
   ["testbed/small/39","small-A01-39"],
   ["testbed/small/40","small-A01-40"],
   ["testbed/small/41","small-A01-41"],
   ["testbed/small/42","small-A01-42"],
   ["testbed/small/43","small-A01-43"],
   ["testbed/small/44","small-A01-44"],
   ["testbed/small/45","small-A01-45"],
   ["testbed/small/46","small-A01-46"],
   ["testbed/small/47","small-A01-47"],
   ["testbed/small/48","small-A01-48"]
]

# define test collections

opt_testbeds = [
   [opt_testsuite_small,opt_testsuite_small_desc]
]

## define architectures to test

opt_testarch = ["CPU","SSE","AVX","AVX2"]


##############################
# DO NOT MODIFY FROM HERE ON #
##############################

colors = {
   "default"  : "",
   "-"        : "\x1b[00m",
   "red"      : "\x1b[31;1m",
   "green"    : "\x1b[32;1m",
   "yellow"   : "\x1b[33;1m",
   "blue"     : "\x1b[34;1m",
   "magenta"  : "\x1b[35;1m",
   "cyan"     : "\x1b[36;1m",
   "bluebg"   : "\x1b[44;1m",
   "yellowbg" : "\x1b[43;2m",
   "cyanbg"   : "\x1b[46;2m",
   "magbg"    : "\x1b[45;2m"
 }

def header():
  sys.stdout.write(" _                   _  _   \n"
                   "| |                 | || |  \n"
                   "| |__  _ __  _ __   | || |_ \n"
                   "| '_ \| '_ \| '_ \  |__   _|\n"
                   "| |_) | |_) | |_) |    | |  \n"
                   "|_.__/| .__/| .__/     |_|  \n"
                   "      | |   | |             \n"
                   "      |_|   |_|             \n"
                   "\n"
                   "bpp 4 testing framework\n\n");
  
def has_colors(stream):
  if not hasattr(stream,"isatty"):
    return False
  if not stream.isatty():
    return False
  try:
    import curses
    curses.setupterm()
    return curses.tigetnum("colors") > 2
  except:
    return False

def ansiprint(color,text,breakline=0):
  if colors[color] and has_colors:
    sys.stdout.write(colors[color] + text + "\x1b[00m")
  else:
    sys.stdout.write(text)
  if breakline:
    sys.stdout.write("\n")

def test_ok():
  ansiprint("green","OK")

def test_fail():
  ansiprint("red","Fail")
  

def difftest(test):
  mcmcfile = test + "/out/mcmc.txt"
  outfile  = test + "/out/out.txt"

  refmcmcfile = test + "/ref/mcmc.txt"

  p = Popen(["diff","-q",mcmcfile,refmcmcfile], stdout=PIPE)
  output = p.communicate()[0]
  return output

def testf(curtest,numtest,t,desc,arch):
  ctl = t + "/data/bpp.ctl"
  cmd = opt_bpp_bin + " --cfile " + ctl + " --arch "  + arch + " 2>tmperr >tmp"

  now = time.strftime("  %H:%M:%S")

  tstart = time.time()

  p1 = Popen(cmd, shell=True)
  os.waitpid(p1.pid,0)

  tend = time.time()
  runtime = tend - tstart

  ansiprint("-", "{:>3}/{:<3} ".format(curtest,numtest) + now)
  ansiprint("cyan", " {:<39} ".format(desc))

  runtime = "%.2f" % runtime
  result = difftest(t)
  ansiprint("cyan", "{:<14} ".format(runtime))
  if result == "":
    test_ok()
  else:
    test_fail()
  print
   
def runtests():
  total = 0;
  for tb in opt_testbeds:
    testlist = tb[0]
    desc = tb[1]
    total += len(testlist) 
  print(" %d tests found" % total)
  print(" %d arch sets" % len(opt_testarch))

  for arch in opt_testarch:
    ansiprint("bluebg", "{:<80}"
               .format(arch.rjust(40+len(arch)/2)), True)
    ansiprint("yellowbg", "{:<7}   {:<8} {:<39} {:<14} Result"
               .format(" ","Start", "Test", "Time [s]"),True)

    current = 0;
    for tb in opt_testbeds:
      testlist = tb[0]
      desc = tb[1]
      for t in testlist:
        test = t[0]
        testdesc = t[1];
        current = current+1
        testf(current,total,test,testdesc,arch)

if __name__ == "__main__":
  
  header()

  runtests()
