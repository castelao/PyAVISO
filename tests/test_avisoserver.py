# -*- coding: utf-8 -*-

import urllib

def test_avisoserver():
    """ Test is AVISO DAP server is still running
    """
    
    urlbase = 'http://opendap.aviso.oceanobs.com/thredds/catalog.html'
    content = urllib.urlopen(urlbase).read()
    assert 'Delayed Time Data' in content
    assert 'Near-Real Time Data' in content
