import os

def create_ff(x, y):
  ff = f'''setParameter[projection=EPSG:32632]
setParameter[fuelsTableFile=./fuels.ff]
setParameter[spatialIncrement=3]
setParameter[propagationModel=Rothermel]
setParameter[minSpeed=0.005]
setParameter[dumpMode=json]
setParameter[caseDirectory=.]
setParameter[ForeFireDataDirectory=.]
setParameter[propagationSpeedAdjustmentFactor=1]
loadData[landscape.nc;2009-07-24T11:37:39Z]
startFire[loc=({x},{y},0);t=0]
step[dt=12000]
print[./*count*-*ISOdate*.json]
print[]
  '''

  dir_path = os.path.dirname(os.path.realpath(__file__))

  c = [str(round(c, 4)).replace('.','_') for c in [x, y]]
  filename = f'{c[0]}__{c[1]}.ff'
  
  output_path = dir_path + '/../examples/aullene/'
  complete_path = output_path + filename

  with open(complete_path, 'w', encoding='utf-8') as f:
    f.write(ff)

  return {'output_path': output_path, 'filename': filename}