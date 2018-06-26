""" Script to read a runs description file and start serial runs on
specified hosts """

import xml.dom.minidom
import os
import stat
import sys
import getopt
import subprocess
import math

class RunsParser:
    " Parse runs definition XML to recover run parameters "

    def __init__(self, runsfilename):
        "Reads XML file, parses it and calls setRuns to store run nodes"

        try:
            filehandle = open(runsfilename, "r")
            try:
                stream = filehandle.read()
            finally:
                filehandle.close()
        except IOError:
            print "File %s doesn't exist\n" % runsfilename
            sys.exit(2)

        runsdoc = xml.dom.minidom.parseString(stream)
        runsnode = runsdoc.documentElement
        self.setRuns(runsnode)

    def setRuns(self, runsnode):
        "stores a dictionary of run names to run XML node objects "

        runelems = runsnode.getElementsByTagName("Run")
        self.runs = {}

        for runelem in runelems:
            nameattr = runelem.attributes["name"]
            runname = nameattr.value
            self.runs[runname] = runelem

    def getRunNames(self):
        " returns a list of run names "

        return self.runs.keys()

    def getType(self, runname):
        """ For the specified run name, return the type of run "SSH", or
        "PBS" """

        runelem = self.runs[runname]

        typeattr = runelem.attributes["type"]

        return typeattr.value


    def getHostsProcs(self, runname):
        """ For the specified run name, return a dictionary of host 
            names and number of processors """

        runelem = self.runs[runname]

        hostelems = runelem.getElementsByTagName("Host")

        if not len(hostelems):
            print "No Hosts found for run %s" % runname
            sys.exit(2)

        hostprocs = {}

        for hostelem in hostelems:

            hostname = hostelem.attributes["name"].value 
            numprocs = hostelem.attributes["numprocs"].value
            hostprocs[hostname] = numprocs

        return hostprocs

    def getQueueProcs(self, runname):
        """ For the specified runname, return a dictionary of queue name
        and the number of processors """

        queuedesc = self.getAttributesOfUniqueElement(runname, "Queue")

        return queuedesc

    def getComm(self, runname, commname):
        " For the specified run name, return the command string "

        runelem = self.runs[runname]

        commelem = runelem.getElementsByTagName(commname)

        if len(commelem) != 1:
            print "More than one %s element in %s!" % (commname, runname)
            sys.exit(2)

        return commelem[0].attributes["command"].value 

    def getAracneOptions(self, runname):
        " For the specified run name, return the command string "

        runelem = self.runs[runname]

        commelem = runelem.getElementsByTagName("aracne")

        if len(commelem) != 1:
            print "More than one %s element in %s!" % ("aracne", runname)
            sys.exit(2)

        return commelem[0].attributes["options"].value 

    def getAracneFile(self, runname):
        runelem = self.runs[runname]
        commelem = runelem.getElementsByTagName("aracne")
        return commelem[0].attributes["filename"].value

    def getCircuitDesc(self, runname):
        " For the specified run, return a dictionary of Circuit attributes "

        runelem = self.runs[runname]

        circuitelem = runelem.getElementsByTagName("Circuit")

        if len(circuitelem) != 1:
            print "More than one Circuit element in %s!" % runname
            sys.exit(2)

        circattrs = circuitelem[0].attributes
        circdesc = {}

        for circkey in circattrs.keys():
            circdesc[circkey] = circattrs[circkey].value


        return circdesc

    def getAttributesOfUniqueElement(self, runname, elemname):
        " For the specified run, return a dictionary of Circuit attributes "

        runelem = self.runs[runname]

        elem = runelem.getElementsByTagName(elemname)

        if len(elem) != 1:
            print "More than one %s element in %s!" % (elemname, \
                                                       runname)
            sys.exit(2)

        attrs = elem[0].attributes
        attrsdict = {}

        for key in attrs.keys():
            attrsdict[key] = attrs[key].value


        return attrsdict



class Runs:
    """ Compile the host information and commands to write run scripts, and
    optionally to dispatch the scripts. """

    def __init__(self, runsdesc, runname):
        " Compile run information based on the runsdesc RunsParser object "

        self.runfilehdr_SSH = "#!/bin/bash"

        self.runfilebody_SSH = """for i in %s;
    do
        %s 

    done	
"""
        self.runfilehdr_PBS = """#!/bin/bash
#PBS -N %s
#PBS -t %s
#PBS -d %s
#PBS -q %s
#PBS -j oe
#PBS -S /bin/bash
"""

        self.runfilebody_PBS = """if (( ${PBS_ARRAYID} == %d ));
then
    for i in %s;
    do
        %s 
    done	
fi
"""
        self.runsdesc = runsdesc
        self.runname = runname

        self.setRunInfo()

        " Get the methods "
        setDispatchInfoMethod = getattr(self, \
                                       "setDispatchInfo_%s" % self.type)

        writeRunScriptMethod = getattr(self, \
                                       "writeRunScript_%s" % self.type)

        dispatchMethod = getattr(self, \
                                 "dispatch_%s" % self.type)


        " Execute "
        setDispatchInfoMethod()
        writeRunScriptMethod()
        dispatchMethod()


    def setRunInfo(self):
        " Set run information "

        self.type = self.runsdesc.getType(self.runname)
        circdesc = self.runsdesc.getAttributesOfUniqueElement(\
                                                    self.runname, \
                                                    "Circuit")

        self.wd = os.path.expanduser(circdesc["dirname"])

        self.circsperproc = int(circdesc["circsperproc"])
        self.startnum = int(circdesc["startnum"])
        self.endnum = int(circdesc["endnum"])
        self.numcircdigits = len(circdesc["startnum"])

        """ prepare aracne and transc commands and state and time file
            names """
        ext = circdesc["extension"]
        stub = circdesc["stubname"]
        self.nicelevel = circdesc["nicelevel"]

        aracnecommdesc = self.runsdesc.getAttributesOfUniqueElement(\
                                                    self.runname, \
                                                    "aracne")

        aracnecomm = os.path.expanduser(aracnecommdesc["command"])
        aracnecomm = aracnecomm + " " + \
                                            aracnecommdesc["options"]
        self.aracnecomm = aracnecomm+" "+"-s "+stub+"_${i}."+ext+" -o "+stub+"_${i}.adj"
        transccommdesc = self.runsdesc.getAttributesOfUniqueElement(\
                                                    self.runname, \
                                                    "Transc")

        transccomm = os.path.expanduser(transccommdesc["command"])
        transccomm = transccomm + " " + \
                                            transccommdesc["options"]
        transccomm = "nice -n "+self.nicelevel+" "+transccomm
        transccomm = transccomm+" "+stub+"_${i}_${j}."+ext
        self.transccomm = transccomm+" &> "+stub+"_${i}_${j}."+ext+".err"

        self.statefile = stub+"_${i}_${j}."+ext+".state"
        self.timesfile = stub+"_${i}_${j}."+ext+".times"

        " Calculate the no of digits for replicate numbers "
        numreplicates = int(circdesc["numreplicates"])
        numrepldigits = len(circdesc["numreplicates"])

        " make the list of replicate numbers w/padding to loop over "
        offset = circdesc.has_key("replicateoffset") and \
                    int(circdesc["replicateoffset"]) or \
                    0
        self.jlooplist = ' '.join([str(offset+j+1).zfill(numrepldigits) \
                                            for j in range(numreplicates)])


        self.startnumlist = range(self.startnum, self.endnum+1, \
                                                self.circsperproc)

    def getEnd(self, start):

        return (start < self.endnum - self.circsperproc) \
                                and start+self.circsperproc \
                                or self.endnum + 1

    def getCircLoopList(self, start):
        " make the list of circuit numbers w/padding to loop over "

        end = self.getEnd(start)
        return ' '.join([str(i).zfill(self.numcircdigits) \
                                        for i in range(start, end)])
   
    def getRunFilename_SSH(self, start):        
        " Get name of SSH run script file "

        end = self.getEnd(start)
        return "run_"+\
                        str(start).zfill(self.numcircdigits)+"-"+\
                        str(end-1).zfill(self.numcircdigits)

    def getRunFilename_PBS(self):        
        " Get name of PBS run script file "

        return "run_"+\
                    str(self.startnum).zfill(self.numcircdigits)+"-"+\
                    str(self.endnum).zfill(self.numcircdigits)


    def writeRunScript_SSH(self):
        " Write a run script to be dispatched by SSH "

        for start in self.startnumlist:
            
            runfiletextlist = [self.runfilehdr_SSH]
            ilooplist = self.getCircLoopList(start)
            " compile data to write to run script "
            runfiletextlist.append(self.runfilebody_SSH % (\
                                        ilooplist,\
                                        self.aracnecomm))

            runfiletext = "\n\n".join(runfiletextlist)

            filename = self.getRunFilename_SSH(start)
            absname = os.path.join(self.wd, filename)
            self.writeFile(absname, runfiletext)

            os.chmod(absname, \
                             stat.S_IRWXU | \
                             stat.S_IRGRP | stat.S_IXGRP | \
                             stat.S_IROTH | stat.S_IXOTH)

    def writeRunScript_PBS(self):
        " Write a run script to be dispatched by SSH "

        arrayjoblist = range(len(self.startnumlist))
        jobnumrange = len(arrayjoblist)-1 and \
                      "-".join([str(arrayjoblist[0]),\
                                str(arrayjoblist[-1])]) or \
                      "0" 
        
        runfiletextlist = [self.runfilehdr_PBS % (\
                                        self.runname,\
                                        jobnumrange,\
                                        self.queuename,\
                                        self.wd)]

        for start in self.startnumlist:

            jobnum = arrayjoblist.pop()
            ilooplist = self.getCircLoopList(start)
            " compile data to write to run script "
            runfiletextlist.append(self.runfilebody_PBS % (\
                                        jobnum,\
                                        ilooplist,\
                                        self.aracnecomm))

           
        runfiletext = "\n\n".join(runfiletextlist)

        filename = self.getRunFilename_PBS()
        absname = os.path.join(self.wd, filename)
        self.writeFile(absname, runfiletext)


    def writeFile(self, filename, data):
        " write out the script file "

        self.writeAracneProbes()

        try:
            filehandle = open(filename, "w")
            try:
                filehandle.write(data)
            finally:
                filehandle.close()
        except IOError:
            print "Couldn't write to file %s\n" % filename
            sys.exit(2)


    def setDispatchInfo_SSH(self):
        """ Check that we have enough proc. Build a list of host names, 
            each host appearing as many times as the number of available 
            procs """

        availhostsprocs = self.runsdesc.getHostsProcs(self.runname)

        " Calculate number of required procs "
        procsreq = math.ceil(\
                    (float(self.endnum)-float(self.startnum)+1.0)/ \
                    float(self.circsperproc))

        availprocs = sum([int(p) for p in availhostsprocs.values()])

        if availprocs < procsreq:
            print "Need %d procs, only %d available!" % \
                                        (procsreq, availprocs)
            sys.exit(2)                                        

        hostlist = []
        for host in availhostsprocs.keys():
            hostlist.extend([host for j in \
                                range(int(availhostsprocs[host]))])

        self.hostlist = hostlist

    def setDispatchInfo_PBS(self):
        """ Check that we have enough procs. Set the queue name. """

        queuedesc = self.runsdesc.getQueueProcs(self.runname)

        " Calculate number of required procs "
        procsreq = math.ceil(\
                    (float(self.endnum)-float(self.startnum)+1.0)/ \
                    float(self.circsperproc))

        availprocs = int(queuedesc["numprocs"])

        if availprocs < procsreq:
            print "Need %d procs, only %d available!" % \
                                        (procsreq, availprocs)
            sys.exit(2)                                        

        
        self.queuename = queuedesc["name"]


    def dispatch_SSH(self):                
        """ compile the command to ssh into remote host and run the
            script. The bash script is launched with nohup and all IO
            redirected so that ssh terminates after launching the
            script in background. """

        for start in self.startnumlist:

            end = self.getEnd(start)
            host = self.hostlist.pop()

            runfilename = self.getRunFilename_SSH(start)
            dispatch_command = "cd "+self.wd+"; nohup ./"+runfilename+\
                               " > /dev/null 2> /dev/null < /dev/null &"

            print "Dispatching circuits "+\
                            str(start).zfill(self.numcircdigits)+"-"+\
                            str(end-1).zfill(self.numcircdigits)+\
                            " on "+host+"..."
            " Execute "
            try:
                ret = subprocess.Popen(["ssh", "-x", \
                                        host, dispatch_command])
            except OSError:
                print "Trouble executing ssh "
                sys.exit(2)

    def dispatch_PBS(self):                
        " Do nothing to dispatch PBS script. That is launched by hand. "

        return

    def writeAracneProbes(self):
        " writes the probes "

        
        filename = self.runsdesc.getAracneFile(self.runname)
        queue = self.runsdesc.getQueueProcs(self.runname)
        procnum = int(queue["numprocs"])

        #generate probe files for aracne
        itemcount = 0
        filenum = 0
        filesize = 0
        filecount = 1

        os.system("cut -f1 {}.txt | awk 'NR>=2' > {}_probes.txt".format(filename, filename))
        with open ("{}_probes.txt".format(filename)) as f:
            probelist = [x.strip() for x in f.readlines()]
            filesize = len(probelist) / procnum
            for i in probelist:
                if not itemcount % filesize:
                    outfile = open(filename+"_"+str(filecount).zfill(8)+".txt", "w")
                    filecount += 1
                filenum += 1
                outfile.write(i+"\n")
                itemcount += 1
            outfile.close()

if __name__ == "__main__":

    runsfilename = ''

    try:	
        (opts, args) = getopt.getopt(sys.argv[1:], "f:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-f":
            runsfilename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    if len(args) != 0:
        print """Need no arguments.\n"""
        sys.exit(2)

    runsdesc = RunsParser(runsfilename)
    runnames = runsdesc.getRunNames()
    for runname in runnames:
        run = Runs(runsdesc, runname)
