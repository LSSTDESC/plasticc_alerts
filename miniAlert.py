# This code includes code from tests directory of the ap_association library
import os
import numpy as np
import pandas as pd
import shutil
import tempfile
import unittest
import datetime
from lsst.ap.association import (PackageAlertsTask,
                                 AssociationTask,
                                 make_dia_source_schema,
                                 make_dia_object_schema)
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.fits as afwFits
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.daf.base as dafBase
from lsst.dax.apdb import Apdb, ApdbConfig
import lsst.geom as geom
import lsst.meas.base.tests
from lsst.utils import getPackageDir
import lsst.utils.tests
import lsst.afw.geom as afwGeom


def makeExposure(flipX=False, flipY=False):
    """Create an exposure and flip the x or y (or both) coordinates.

    Returns bounding boxes that are right or left handed around the bounding
    polygon.

    Parameters
    ----------
    flipX : `bool`
        Flip the x coordinate in the WCS.
    flipY : `bool`
        Flip the y coordinate in the WCS.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        Exposure with a valid bounding box and wcs.
    """
    metadata = dafBase.PropertySet()

    metadata.set("SIMPLE", "T")
    metadata.set("BITPIX", -32)
    metadata.set("NAXIS", 2)
    metadata.set("NAXIS1", 1024)
    metadata.set("NAXIS2", 1153)
    metadata.set("RADECSYS", 'FK5')
    metadata.set("EQUINOX", 2000.)

    metadata.setDouble("CRVAL1", 215.604025685476)
    metadata.setDouble("CRVAL2", 53.1595451514076)
    metadata.setDouble("CRPIX1", 1109.99981456774)
    metadata.setDouble("CRPIX2", 560.018167811613)
    metadata.set("CTYPE1", 'RA---SIN')
    metadata.set("CTYPE2", 'DEC--SIN')

    xFlip = 1
    if flipX:
        xFlip = -1
    yFlip = 1
    if flipY:
        yFlip = -1
    metadata.setDouble("CD1_1", xFlip * 5.10808596133527E-05)
    metadata.setDouble("CD1_2", yFlip * 1.85579539217196E-07)
    metadata.setDouble("CD2_2", yFlip * -5.10281493481982E-05)
    metadata.setDouble("CD2_1", xFlip * -8.27440751733828E-07)

    wcs = afwGeom.makeSkyWcs(metadata)
    exposure = afwImage.makeExposure(
        afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))), wcs)
    detector = DetectorWrapper(id=23, bbox=exposure.getBBox()).detector
    visit = afwImage.VisitInfo(
        exposureId=1234,
        exposureTime=200.,
        date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                              dafBase.DateTime.Timescale.TAI))
    exposure.setDetector(detector)
    exposure.getInfo().setVisitInfo(visit)
    exposure.setFilter(afwImage.Filter('g'))

    return exposure


def makeDiaObjects(nObjects, exposure, pixelator):
    """Make a test set of DiaObjects.

    Parameters
    ----------
    nObjects : `int`
        Number of objects to create.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create objects over.
    pixelator : `lsst.sphgeom.HtmPixelization`
        Object to compute spatial indicies from.

    Returns
    -------
    diaObjects : `pandas.DataFrame`
        DiaObjects generated across the exposure.
    """
    bbox = geom.Box2D(exposure.getBBox())
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nObjects)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nObjects)

    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    wcs = exposure.getWcs()

    data = []
    for idx, (x, y) in enumerate(zip(rand_x, rand_y)):
        coord = wcs.pixelToSky(x, y)
        htmIdx = pixelator.index(coord.getVector())
        newObject = {"ra": coord.getRa().asDegrees(),
                     "decl": coord.getRa().asDegrees(),
                     "radecTai": midPointTaiMJD,
                     "diaObjectId": idx,
                     "pixelId": htmIdx,
                     "pmParallaxNdata": 0,
                     "nearbyObj1": 0,
                     "nearbyObj2": 0,
                     "nearbyObj3": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            newObject["%sPSFluxNdata" % f] = 0
        data.append(newObject)

    return pd.DataFrame(data=data)


def makeDiaSources(nSources, diaObjectIds, exposure, pixelator):
    """Make a test set of DiaSources.

    Parameters
    ----------
    nSources : `int`
        Number of sources to create.
    diaObjectIds : `numpy.ndarray`
        Integer Ids of diaobjects to "associate" with the DiaSources.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create sources over.
    pixelator : `lsst.sphgeom.HtmPixelization`
        Object to compute spatial indicies from.

    Returns
    -------
    diaSources : `pandas.DataFrame`
        DiaSources generated across the exposure.
    """
    bbox = geom.Box2D(exposure.getBBox())
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nSources)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nSources)
    rand_ids = diaObjectIds[np.random.randint(len(diaObjectIds), size=nSources)]

    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    wcs = exposure.getWcs()

    data = []
    for idx, (x, y, objId) in enumerate(zip(rand_x, rand_y, rand_ids)):
        coord = wcs.pixelToSky(x, y)
        htmIdx = pixelator.index(coord.getVector())
        data.append({"ra": coord.getRa().asDegrees(),
                     "decl": coord.getRa().asDegrees(),
                     "diaObjectId": objId,
                     "diaSourceId": idx,
                     "pixelId": htmIdx,
                     "midPointTai": midPointTaiMJD})

    return pd.DataFrame(data=data)

hi=makeExposure()