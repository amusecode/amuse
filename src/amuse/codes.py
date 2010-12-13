
import sys

_CODES = [
'athena', 'capreole', 'cachedse',     'gadget2',      'mesa',         'octgrav',      'twobody',
'capreole',     'hermite0',     'mocassin',     'phiGRAPE',
'athena',       'evtwin',       'hop',          'seba',
'bhtree',       'evtwin2sse',   'interface',    'smallN',
'bse',          'fi',           'mercury',      'sse',
]

__all__ = []

def _import_modules():

    for x in _CODES:
        modulename = 'amuse.legacy.' + x + '.interface'
        try:
            __import__(modulename)
            globals()[x] = sys.modules[modulename]
            __all__.append(x)
        except ImportError as ex:
            modulename = 'amuse.legacy.' + x + '.' + x
            try:
                __import__(modulename) 
                globals()[x] = sys.modules[modulename]
                __all__.append(x)
            except ImportError as ex:
                pass

_import_modules()
    