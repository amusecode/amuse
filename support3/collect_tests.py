from nose.plugins import Plugin
import collections


class CollectTests(Plugin):
    """collect test cases and output json"""
    
    name = 'collecttests'
    score = 1100

    useDetails = False
    config = None
    loaderClass = None
    
    def options(self, parser, env):
        """Register commandline options
"""
        Plugin.options(self, parser, env)

    def configure(self, options, conf):
        if not self.can_configure:
            return
        Plugin.configure(self, options, conf)
        
        self.config = conf

    #copied from multiprocess plugin
    def prepareTestLoader(self, loader):
        """Remember loader class so MultiProcessTestRunner can instantiate
the right loader.
"""
        self.loaderClass = loader.__class__

        #monkey patch SuiteFactory.__call__ function to count how many tests there are
        sC = loader.suiteClass
        plugin = self
        def __call__(tests=(), **kw):
            is_suite = isinstance(tests, TestSuite)
            if isinstance(tests, collections.Callable) and not is_suite:
                tests = list(tests())
            elif is_suite:
                tests = [tests]
            suite = sC(tests, **kw)
            if not is_suite:
                if not getattr(tests, '__len__', None):
                    lt = tests.tests
                else:
                    lt = tests
                #preload all lazy load suite, so we can count totalTests properly
                if self.preload:
                    for t in lt:
                        plugin.preloadLazySuite(t)
                plugin.totalTests += sum(1 for i in lt if not isinstance(i, TestSuite))
                #plugin.totalTests += sum(i.countTestCases() for i in lt)
                #import pdb;pdb.set_trace()
                #print "totalTests", plugin.totalTests, '*'*60
            return suite
        #assigning to suiteClass.__call__ won't make it work: calling suiteClass() won't trigger the __call__ attribute
        loader.suiteClass = __call__
    def preloadLazySuite(self, t):
        if isinstance(t, LazySuite):
            tests = [i for i in t]
            t._tests = tests
    #so we stick with prepareTestResult for now
    def prepareTestRunner(self, runner):
        #replace _makeResult in the default nose TestRunner to return
        #our implementation SubunitTestResult

        if not hasattr(runner, "_makeResult"):
            raise Exception('''runner does not have _makeResult method,\
don't know how to attach to it.''')
        
        #this plugin runs before multiprocess plugin, and this function
        #return a runner, so it will prevent multiprocess.prepareTestRunner
        #from executing. so if multiprocess is enabled, we need to create
        #MultiProcessTestRunner
        if self.multiprocess_workers:
            #with MultiProcessTestRunner.run, beforeTest is called too
            #early, so it will never print any progress info
            #instead, we have to monkey patch nextBatch which is
            #called by MultiProcessTestRunner.run to collect all tests
            from nose.plugins.multiprocess import MultiProcessTestRunner
            class SubunitMultiProcessTestRunner(MultiProcessTestRunner):
                _top = 1
                def nextBatch(self, test):
                    top = self._top
                    if top:
                        self._top = 0
                    for batch in MultiProcessTestRunner.nextBatch(self, test):
                        yield batch
                    if top:
                        plugin.printProgress()

            runner = SubunitMultiProcessTestRunner(stream=runner.stream,
                                      verbosity=self.config.verbosity,
                                      config=self.config,
                                      loaderClass=self.loaderClass)

        runner.useDetails = self.useDetails
        plugin = self
        def _makeResult(self):
            result = SubunitTestResult(self.stream, self.descriptions,
                self.config, useDetails=self.useDetails)
            result.addTime()
            #save the result so it can be used in beforeTest below
            plugin.result = result
            #plugin.printProgress()
            #if plugin.totalTests:
            # result.progress(plugin.totalTests, PROGRESS_CUR)
            # result.addTime()
            return result
        runner._makeResult = instancemethod(_makeResult,
          runner, runner.__class__)

        return runner
    totalTests = 0
    def printProgress(self):
        if self.totalTests:
            self.result.progress(self.totalTests, PROGRESS_CUR)
            self.totalTests = 0
    def beforeTest(self, *args):
        self.printProgress()
