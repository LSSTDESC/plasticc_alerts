import io
from ligo.skymap.io import read_sky_map
from lsst.alert.packet import Schema
import numpy as np

import ligo.skymap.plot
import matplotlib.pyplot as plt

skymap_filename = 'plasticc_schema/sample_data/S190814bv-bayestar.multiorder.fits'

skymap, header = read_sky_map(skymap_filename)
with open(skymap_filename, 'rb') as f:
    skymap_bytes = f.read()

schema = Schema.from_file('plasticc_schema/lsst.v4_1.lvkAlertContent.avsc')

mock_alert_content = dict(
    supereventId='S123456',
    gpstime=header['gps_time'],
    skymapFilename='bayestar.fits.gz',
    skymapHealpix=skymap_bytes
)
serialized_alert = schema.serialize(mock_alert_content)

# sanity check
deserialized_alert = schema.deserialize(serialized_alert)
hpx, header = read_sky_map(io.BytesIO(deserialized_alert['skymapHealpix']))
import healpy
healpy.mollview(hpx)
plt.show()
