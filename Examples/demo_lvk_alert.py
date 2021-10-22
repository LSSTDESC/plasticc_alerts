from ligo.skymap.io import read_sky_map
from lsst.alert.packet import Schema
import numpy as np

import ligo.skymap.plot
import matplotlib.pyplot as plt

skymap, header = read_sky_map('plasticc_schema/sample_data/S190814bv-bayestar.multiorder.fits')
schema = Schema.from_file('plasticc_schema/lsst.v4_1.lvkAlertContent.avsc')

mock_alert_content = dict(
    supereventId='S123456',
    gpstime=header['gps_time'],
    skymapFilename='bayestar.fits.gz',
    skymapHealpix=skymap
)
serialized_alert = schema.serialize(mock_alert_content)

# sanity check
deserialized_alert = schema.deserialize(serialized_alert)
plt.figure()
ax = plt.axes(projection='astro globe')
ax.imshow_hpx(np.array(deserialized_alert['skymapHealpix']))
plt.show()
