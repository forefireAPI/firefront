import os
from datetime import date
from .ff2geojson import *

class Forefire:

  version = '2.0.0'

  def __init__(self):
    self.ff = '''setParameter[dumpMode=geojson]
setParameter[caseDirectory=.]
setParameter[ForeFireDataDirectory=.]
'''
    self.year = date.today().year
    self.month = date.today().month
    self.day = date.today().day

  def setProjection(self, proj='EPSG:32632'):
    self.ff += f'setParameter[projection={proj}]\n'

  def setFuels(self, fuelsTableFile='./fuels.ff'):
    self.ff += f'setParameter[fuelsTableFile={fuelsTableFile}]\n'

  def setPropagationModel(self, propagationModel='Rothermel'):
    self.ff += f'setParameter[propagationModel={propagationModel}]\n'

  def setDate(self):
    self.ff += f'setParameter[year={self.year}]\n'
    self.ff += f'setParameter[month={self.month}]\n'
    self.ff += f'setParameter[day={self.day}]\n'

  def setParameter(self, parameter, value):
    self.ff += f'setParameter[{parameter}={value}]\n'

  def loadData(self, nc='landscape.nc', isoDate='2009-07-24T11:37:39Z'):
    self.ff += f'loadData[{nc};{isoDate}]\n'

  def setFireDomain(self, sw, ne):
    self.ff += f'FireDomain[sw={sw};ne={ne};t=0.]\n'

  def setFirefront(self, coords_list, vel_list):
    if len(coords_list) == len(vel_list) > 0:
      self.ff += 'FireFront[t=0.]\n'
      for i in range(len(coords_list)):
        self.ff += f'FireNode[loc={coords_list[i]};vel={vel_list[i]};t=0.]\n'

  def startFire(self, lon=499073.45383159, lat=4619272.9498144, t=0):
    self.ff += f'startFire[loc=({lon},{lat},0);t={t}]\n'

  def goTo(self, t):
    self.ff += f'goTo[t={t}]\n'

  def step(self, dt):
    self.ff += f'step[dt={dt}s]\n'

  def printOutput(self):
    self.ff += 'print[./*count*-*ISOdate*.ffgeojson]\nprint[]'

  def saveFf(self, path):
    with open(path, 'w', encoding='utf-8') as f:
      f.write(self.ff)
  
  def configBasicFf(self, lon, lat, t=12000):
    self.setProjection()
    self.setFuels()
    self.setPropagationModel()
    self.setDate()
    self.loadData()
    self.startFire(lon, lat)
    self.step(t)
    self.printOutput()

  def convert_to_geojson(self, path):
    ffjson2geojson(path)
    return True