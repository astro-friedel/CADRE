import calculations

def writeLog(project,source,flux,uncert,frequency,SB,obsDate,primary):
    #fl = open("/export/scr0/pipeline/global-scratch/input/test/flux.html",'r')
    #lines = fl.readlines()
    #fl.close()

    date = calculations.numericToConventional(obsDate)

    fl = open("/export/drobo3/pipeline/test/flux.html",'a')
#    print len(lines)
#    print lines
    #for i in range(0,len(lines)-2) :
    #    fl.write(lines[i])
#        print lines[i]," ",i
    #format is: source, date, project, sideband, primary flux cal, frequency, flux (uncertainty)
    outline = "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.2f</td><td>%.2f (%.3f)</td></tr>\n" % (source,date,project,SB,primary,frequency,flux,uncert)
    fl.write(outline)
#    print outline
#    print lines[-2]," -2\n"
#    print lines[-1]," -1\n"
    #fl.write(lines[-2])
    #fl.write(lines[-1])
    fl.close()
