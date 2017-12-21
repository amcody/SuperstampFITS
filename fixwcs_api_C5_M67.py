import os
import myclient
import time

files = open('filenames_M67.txt', 'r')
for line in files:
    line = line.strip()
    end = line.find('.fits')
    start = line.find('M67_')
    print line

    cli = myclient.MyClient()

### Put your API key below in places of xxxxxxx
    cli.login('xxxxxxxxxxxxxx')

    #Set some solving options
    opts=dict()
    opts['allow_commercial_use'] = 'n'
    opts['allow_modifications'] = 'n'
    opts['publicly_visible'] =  'n'
    opts['scale_units'] = "arcsecperpix"
    opts['scale_lower'] = 3.9
    opts['scale_upper'] = 4.1
    opts['tweak_order'] = 1

    #Upload an image
    cli.upload(line, **opts)

    #Save state in case of problems. This uses python shelves underneath
    #so it's not portable between machines. To restore a client in
    #a new session, use cli.load()
    cli.save('cli.sav')

    data = cli.getFilenames()
    job = data[0][0]

# See if job is done; if not, sleep a bit:
    status = cli.check_status_single_job(job)
    check = (status['processing_finished'] != 'None') and ('success' in (status['jobStatus']).values())
    failed = ('failure' in (status['jobStatus']).values())

    while (check == False and failed == False):
       time.sleep(5)                    
       print 'Not done yet'
       try:
         status = cli.check_status_single_job(job)
         check = (status['processing_finished'] != 'None') and ('success' in (status['jobStatus']).values())
         failed = ('failure' in (status['jobStatus']).values())
         print check
       except:
         print 'Error- trying again'

    if (failed != True):
      newfile = line[start:end]+'_wcs.fits'
      print newfile
      cli.getNewImage(job,newfile)

files.close()
