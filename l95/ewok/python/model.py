import os
import ewok

__all__ = ["fcfile", "fbfile"]


def fcfile(fcout, step):
    fc = {}
    fc['date'] = fcout['date']
    keys = [fcout['exp'], fcout['type'], fcout['date'], ewok.jediformat(step)]
    fname = '.'.join(keys)
    fc['filename'] = os.path.join(fcout['datadir'], fname)
    return fc

#def obs_file(conf):
#    obsfile = conf['obsdatain']
#    return obsfile

def fbfile(conf):
    fbfile = conf['obsdataout']
    return fbfile

