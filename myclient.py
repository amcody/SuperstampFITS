import client
import urllib2
import simplejson
import re
import shelve
import pyfits

## Put your API key here:
apiKey="xxxxxxxxxxxxxxxxx"



class MyClient(client.Client):
    def __init__(self, apiurl=None):
        """A interface to the Astrometry.net SOAP protocol

        Extends the client class provided by the project
        that can be found at
        http://astrometry.net/svn/trunk/src/astrometry/net/client

        The python files at that URL need to be checked out an
        in your PYTHONPATH for this client to work.

        See example() in this module for example code
        """

        if apiurl is None:
            apiurl = client.Client.default_url

        client.Client.__init__(self, apiurl)
        self.submittedJobs = dict()
        self.jobStatus = dict()


    def upload(self, fn, **kwargs):
        """Upload a file to astrometry.net server and solve it

        Input:
        fn  (string)    Name of file to upload

        Optional arguments.
        Optional arguments aren't documented anywhere(?) See example()
        in this module for some known ones

        Returns:
        Id number of this submission. Use this id to track a submission
        """
        ret = client.Client.upload(self, fn, **kwargs)

        subId = ret['subid']
        #subId = 1
        self.submittedJobs[subId] = ret
        self.submittedJobs[subId]['filename'] = fn
        self.submittedJobs[subId]['opts'] = kwargs
        del self.submittedJobs[subId]['status']

        return subId

    def save(self, filename):
        """Save a client instance, so you can restore your
        job ids later

        Input:
        filename    (String) Name of file to save session to.

        Notes:
        This method uses the python shelf module, so sessions are not
        portable between machines.
        """

        shelf = shelve.open(filename)
        shelf['submittedJobs'] = self.submittedJobs
        shelf['jobStatus'] = self.jobStatus
        shelf['session'] = self.session
        shelf['apiurl'] = self.apiurl
        shelf.close()


    def load(self, filename):
        """Restore the state of a client from file"""
        shelf = shelve.open(filename)
        self.submittedJobs = shelf['submittedJobs']
        self.jobStatus     = shelf['jobStatus']
        self.session       = shelf['session']
        self.apiurl        = shelf['apiurl']


    def check_status(self):
        """Check status of all jobs associated with this client"""
        for subId in self.submittedJobs.keys():
            self.check_status_single_job(subId)

    def check_status_single_job(self, subId):
        """Check status of a single job

        Input:
        subId (int) id of job as returned, eg by upload()

        Returns
        A dictionary of parameters. The most interesting is
        jobStatus, which will be "success", "failed", or "solving"
        """
        url = self.apiurl + "submissions/%i" %(subId)

        request = urllib2.Request(url)
        print url
        print request
        f = urllib2.urlopen(request)
        txt = f.read()
        result = simplejson.loads(txt)

        for k in result.keys():
            self.submittedJobs[subId][k] = result[k]

        #Now check their status for the various job ids for this submission id
        if "jobStatus" not in self.submittedJobs[subId].keys():
            self.submittedJobs[subId]['jobStatus'] = dict()

        for j in result['jobs']:
            url = self.apiurl + "jobs/%i" %(j)
            request = urllib2.Request(url)
            f = urllib2.urlopen(request)
            sResult = simplejson.loads(f.read())

            self.submittedJobs[subId]['jobStatus'][j] = sResult['status']
            print "Submission: %i Job: %i Status: %s" %(subId, j, sResult['status'])

        return self.submittedJobs[subId]

    def getFilenames(self):
        """Get the names of the files submitted by the client"""
        return self.getValues(['filename'])


    def getValues(self, valueList):
        """Extract a list of properties for each file submitted by the
        client
        """

        if not isinstance(valueList, list):
            valueList = [valueList]

        out = []
        for k in self.submittedJobs.keys():
            row = [k]
            for v in valueList:
                row.append( self.submittedJobs[k][v])
            out.append(row)
        return out


    def getWcs(self, subId, outFilename=None):
        """Download the WCS and return as a pyfits header object.


        Inputs:
        subId
        outFilename (string) If not None, save the header to this file

        Returns:
        A pyfits.Header() object
        """

        jobId = self.submittedJobs[subId]['jobs'][0]
        url = re.sub("api", "wcs_file/%i" %(jobId), self.apiurl)
        #return url
        f = urllib2.urlopen(url)
        text = f.read()

        hdr= pyfits.Header().fromstring(text)

        if outFilename is not None:
            wf = open(outFilename, "wb")
            wf.write(text)
            wf.close()
        return hdr


    def getNewImage(self, subId, outFilename):
        """Download the solved image with WCS keywords in header

        Inputs:
        subId
        outFilename (string) Name of file to save image as

        Returns:
        A pyfits.Header() object
        """
        jobId = self.submittedJobs[subId]['jobs'][0]
        url = re.sub("api", "new_fits_file/%i" %(jobId), self.apiurl)
        f = urllib2.urlopen(url)

        w = open(outFilename, "wb")
        w.write(f.read())
        w.close()
        return True





def example():
    #Create a new client
    cli = myclient.MyClient()

    #Login. Create an API key at http://nova.astrometry.net/api_help
    cli.login(apiKey)

    #Set some solving options
    opts=dict()
    opts['allow_commercial_use'] = 'n'
    opts['allow_modifications'] = 'n'
    opts['publicly_visible'] =  'n'
    opts['scale_units'] = "arcsecperpix"
    opts['scale_lower'] = 3.9
    opts['scale_upper'] = 4.1

    #Upload an image
    cli.upload("file.fits", **opts)

    #Save state in case of problems. This uses python shelves underneath
    #so it's not portable between machines. To restore a client in
    #a new session, use cli.load()
    cli.save('cli.sav')

    #Get a list of identifies for submitted jobs
    data = cli.getFilenames()
    key = data[0][0]

    #Get status
    cli.check_status_single_job(key)
    #or
    print cli.getValues(['filename', 'jobStatus'])

    #Download a WCS
    cli.getWcs(key)

    #Download the new image
    cli.getNewImage(job,'new.fits')
