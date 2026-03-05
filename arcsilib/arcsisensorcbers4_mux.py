"""
Module that contains the ARCSICBERS4MUXSensor class.
"""

import collections
import datetime
import json
import math
import os
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy
import Py6S
import rasterio
import rsgislib
from osgeo import gdal, osr

import arcsilib.arcsiutils

from .arcsiexception import ARCSIException
from .arcsisensor import ARCSIAbstractSensor


class ARCSICBERS4MUXSpectralBandObj(object):
    """
    This is a class to store the information associated with a CBERS-4/MUX image band.
    """

    def __init__(
        self,
        phyBandName=None,
        bandID=None,
        imgRes=0.0,
        wvLenMin=0.0,
        wvLenMax=0.0,
        wvLenCen=0.0,
        respFuncStep=1.0,
        respFunc=None,
    ):
        if respFunc is None:
            respFunc = list()
        self.phyBandName = phyBandName
        self.bandID = bandID
        self.imgRes = imgRes
        self.wvLenMin = wvLenMin
        self.wvLenMax = wvLenMax
        self.wvLenCen = wvLenCen
        self.respFuncStep = respFuncStep
        self.respFunc = respFunc
        self.respFunc6S = None
        self.wvLenMin6S = 0.0
        self.wvLenMax6S = 0.0


class ARCSICBERS4MUXSensor(ARCSIAbstractSensor):
    """
    A class which represents the landsat 8 sensor to read
    header parameters and apply data processing operations.
    """

    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "CBERS4_MUX"
        self.collection_num = 0

        self.band1File = ""
        self.band2File = ""
        self.band3File = ""
        self.band4File = ""
        self.row = 0
        self.path = 0

        self.b1RadMulti = 0
        self.b1CalMax = 0
        self.b2RadMulti = 0
        self.b2CalMax = 0
        self.b3RadMulti = 0
        self.b3CalMax = 0
        self.b4RadMulti = 0
        self.b4CalMax = 0

        self.b1RadAdd = 0.0
        self.b1MaxRad = 0.0
        self.b2RadAdd = 0.0
        self.b2MaxRad = 0.0
        self.b3RadAdd = 0.0
        self.b3MaxRad = 0.0
        self.b4RadAdd = 0.0
        self.b4MaxRad = 0.0

        self.b1CalMin = 0.0
        self.b1CalMax = 0.0
        self.b2CalMin = 0.0
        self.b2CalMax = 0.0
        self.b3CalMin = 0.0
        self.b3CalMax = 0.0
        self.b4CalMin = 0.0
        self.b4CalMax = 0.0

        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizePan = 0.0
        self.gridCellSizeRefl = 0.0
        self.gridCellSizeTherm = 0.0

        self.specBandInfo = dict()

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the CBERS-4/MUX .xml header files
        """
        try:
            inputHeaderPath = Path(inputHeader)
            base_dir = inputHeaderPath.parent
            prefix = inputHeaderPath.stem.rsplit("_", 1)[0]

            self.band5MetaFile = base_dir / f"{prefix}_BAND5.xml"
            self.band6MetaFile = base_dir / f"{prefix}_BAND6.xml"
            self.band7MetaFile = base_dir / f"{prefix}_BAND7.xml"
            self.band8MetaFile = base_dir / f"{prefix}_BAND8.xml"
            self.band5File = base_dir / f"{prefix}_BAND5.tif"
            self.band6File = base_dir / f"{prefix}_BAND6.tif"
            self.band7File = base_dir / f"{prefix}_BAND7.tif"
            self.band8File = base_dir / f"{prefix}_BAND8.tif"

            inputHeader = os.path.abspath(inputHeader)
            self.headerFileName = os.path.split(inputHeader)[1]

            tree = ET.parse(inputHeader)
            root = tree.getroot()
            namespaces = {"ns": "http://www.gisplan.com.br/xmlsat"}

            # Get the sensor info.
            satellite = root.find("ns:satellite", namespaces)
            if satellite is not None:
                name = satellite.find("ns:name", namespaces)
                number = satellite.find("ns:number", namespaces)
                instrument = satellite.find("ns:instrument", namespaces)
                if name is not None and number is not None and instrument is not None:
                    name_text = name.text
                    number_text = number.text
                    instrument_text = instrument.text
                    instrument_channel = instrument.get("channel")
                    if (
                        name_text == "CBERS"
                        and number_text == "4"
                        and instrument_text == "MUX"
                        and instrument_channel == "1"
                    ):
                        self.sensor = "CB4MUX"
            else:
                raise ARCSIException(
                    "Do no recognise the spacecraft and sensor or combination."
                )

            # Get row/path
            path_elem = root.find("ns:image/ns:path", namespaces)
            row_elem = root.find("ns:image/ns:row", namespaces)
            if path_elem is not None and path_elem.text:
                self.path = path_elem.text
            if row_elem is not None and row_elem.text:
                self.row = row_elem.text

            # Get date and time of the acquisition
            viewing_center_elem = root.find("ns:viewing/ns:center", namespaces)
            if viewing_center_elem is not None and viewing_center_elem.text:
                acquisitionTime = viewing_center_elem.text
                # Separate date and time
                data_part, time_part = acquisitionTime.split("T")
                # Separate years, months, days
                acData = data_part.split("-")
                # Separate hours, minutes, seconds (if exists)
                if "." in time_part:
                    time_without_ms, milliseconds = time_part.split(".")
                else:
                    time_without_ms = time_part
                    milliseconds = "0"
                # Separate hours, minutes, seconds
                acTime = time_without_ms.split(":")
                # Separate integer seconds (part before decimal, if exists)
                secsTime = acTime[2].split(".")[0] if "." in acTime[2] else acTime[2]
                # Criate datetime
                self.acquisitionTime = datetime.datetime(
                    int(acData[0]),  # year
                    int(acData[1]),  # month
                    int(acData[2]),  # day
                    int(acTime[0]),  # hour
                    int(acTime[1]),  # minute
                    int(secsTime),  # second
                    int(
                        milliseconds.ljust(6, "0")[:6]
                    ),  # microsecond (fill with zeros)
                )

            # Tier (called level in CBERS)
            tier_elem = root.find("ns:image/ns:level", namespaces)
            if tier_elem is not None and tier_elem.text:
                self.tier = tier_elem.text

            # Orbit Direction
            orbitDirection_elem = root.find("ns:image/ns:orbitDirection", namespaces)
            if orbitDirection_elem is not None and orbitDirection_elem.text:
                self.orbitDirection = orbitDirection_elem.text

            # Absolute Calibration Coefficient (DN -> Rad) #TODO
            coef_elem = root.find(".//ns:image/ns:absoluteCalibrationCoefficient", namespaces)
            coeffs = {}
            if coef_elem is not None and coef_elem.text:
                for band in coef_elem.findall("ns:band", namespaces):
                    name = band.attrib["name"]
                    value = float(band.text)
                    coeffs[name] = value
            self.b5RadMulti = rsgislib.tools.utils.str_to_float(coeffs["5"])
            self.b6RadMulti = rsgislib.tools.utils.str_to_float(coeffs["6"])
            self.b7RadMulti = rsgislib.tools.utils.str_to_float(coeffs["7"])
            self.b8RadMulti = rsgislib.tools.utils.str_to_float(coeffs["8"])

            self.b5RadAdd = 0  # CB4MUX only has multiplicative coefficient
            self.b6RadAdd = 0
            self.b7RadAdd = 0
            self.b8RadAdd = 0

            # self.b1ReflMulti = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_MULT_BAND_1"]
            # )
            # self.b2ReflMulti = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_MULT_BAND_2"]
            # )
            # self.b3ReflMulti = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_MULT_BAND_3"]
            # )
            # self.b4ReflMulti = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_MULT_BAND_4"]
            # )

            # self.b1ReflAdd = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_ADD_BAND_1"]
            # )
            # self.b2ReflAdd = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_ADD_BAND_2"]
            # )
            # self.b3ReflAdd = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_ADD_BAND_3"]
            # )
            # self.b4ReflAdd = rsgislib.tools.utils.str_to_float(
            #     headerParams["REFLECTANCE_ADD_BAND_4"]
            # )

            # self.b1CalMin = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MIN_BAND_1"]
            # )
            # self.b1CalMax = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MAX_BAND_1"]
            # )
            # self.b2CalMin = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MIN_BAND_2"]
            # )
            # self.b2CalMax = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MAX_BAND_2"]
            # )
            # self.b3CalMin = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MIN_BAND_3"]
            # )
            # self.b3CalMax = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MAX_BAND_3"]
            # )
            # self.b4CalMin = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MIN_BAND_4"]
            # )
            # self.b4CalMax = rsgislib.tools.utils.str_to_float(
            #     headerParams["QUANTIZE_CAL_MAX_BAND_4"]
            # )

            # if "CLOUD_COVER" in headerParams:
            #     self.cloudCover = rsgislib.tools.utils.str_to_float(
            #         headerParams["CLOUD_COVER"], 0.0
            #     )
            # if "EARTH_SUN_DISTANCE" in headerParams:
            #     self.earthSunDistance = rsgislib.tools.utils.str_to_float(
            #         headerParams["EARTH_SUN_DISTANCE"], 0.0
            #     )

            # Get Projection
            lon_elem = float(root.find(".//ns:originLongitude", namespaces).text)
            datum = root.find(".//ns:datumName", namespaces).text.upper()
            utmZone = math.floor((lon_elem + 180) / 6) + 1
            lat_elem = root.find(".//ns:boundingBox/ns:UL/ns:latitude", namespaces)
            if lat_elem is not None:
                lat = float(lat_elem.text)
            else:
                raise ARCSIException(f"Do no recognise the latitude: {lat}.")
            utmHem = "s" if lat < 0 else "n"

            self.projNameStr = ""
            if utmZone == 0:
                raise ARCSIException(
                    "CBERS-4/MUX Image is not projected with UTM - contact support as header not currently supported."
                )
            else:
                utmZoneStr = str(utmZone).replace("-", "")
                self.projNameStr = "utm" + utmZoneStr + utmHem

            # Get EPSG Code.
            epsg = None
            if datum in ("WGS84", "WGS 84"):
                # UTM WGS84: 326xx (North), 327xx (South)
                epsg = 32700 + utmZone if utmHem == "s" else 32600 + utmZone
            elif datum in ("SIRGAS2000", "SIRGAS 2000"):
                # UTM SIRGAS 2000 (Brazil)
                # South: 3197x (e.g.: zone 23S → 31983)
                epsg = 31960 + utmZone  # valid for all Brazil
            else:
                raise ARCSIException(f"Do no recognise datum: {datum}.")
            self.epsgCode = epsg
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(self.epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()

            specBandObj = ARCSICBERS4MUXSpectralBandObj()

            # BAND5
            specBandObj.phyBandName = "BAND5"
            specBandObj.bandID = "0"
            specBandObj.imgRes = 20
            specBandObj.wvLenMin = 360  # 450
            specBandObj.wvLenMax = 1000  # 520
            specBandObj.wvLenCen = 485.0
            specBandObj.respFuncStep = 1
            # Values imported from TerraLib
            specBandObj.respFunc = [
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000088,
                0.000177,
                0.000088,
                0.000000,
                0.000000,
                0.000000,
                0.000595,
                0.001190,
                0.000595,
                0.000000,
                0.000072,
                0.000143,
                0.000140,
                0.000136,
                0.000068,
                0.000000,
                0.000132,
                0.000265,
                0.000132,
                0.000000,
                0.000122,
                0.000244,
                0.000239,
                0.000234,
                0.000406,
                0.000578,
                0.000289,
                0.000000,
                0.000324,
                0.000648,
                0.000324,
                0.000000,
                0.000347,
                0.000694,
                0.000685,
                0.000676,
                0.000571,
                0.000466,
                0.002548,
                0.004631,
                0.007909,
                0.011188,
                0.017990,
                0.024792,
                0.035768,
                0.046744,
                0.075321,
                0.103898,
                0.156567,
                0.209236,
                0.270533,
                0.331829,
                0.377406,
                0.422984,
                0.456530,
                0.490076,
                0.523851,
                0.557626,
                0.578408,
                0.599189,
                0.602330,
                0.605470,
                0.606610,
                0.607750,
                0.617594,
                0.627439,
                0.638078,
                0.648717,
                0.657317,
                0.665917,
                0.671048,
                0.676180,
                0.686741,
                0.697303,
                0.703910,
                0.710518,
                0.729951,
                0.749384,
                0.748872,
                0.748361,
                0.759258,
                0.770155,
                0.779815,
                0.789476,
                0.801857,
                0.814237,
                0.826101,
                0.837964,
                0.849017,
                0.860071,
                0.871003,
                0.881934,
                0.897109,
                0.912283,
                0.918861,
                0.925439,
                0.941151,
                0.956862,
                0.969419,
                0.981976,
                0.990988,
                1.0,
                0.994021,
                0.988042,
                0.987170,
                0.986297,
                0.991499,
                0.996701,
                0.983983,
                0.971265,
                0.935965,
                0.900665,
                0.852920,
                0.805175,
                0.776565,
                0.747954,
                0.686618,
                0.625282,
                0.491865,
                0.358447,
                0.269545,
                0.180643,
                0.140602,
                0.100561,
                0.081846,
                0.063132,
                0.052699,
                0.042266,
                0.033774,
                0.025282,
                0.020260,
                0.015238,
                0.012193,
                0.009148,
                0.007792,
                0.006436,
                0.005440,
                0.004443,
                0.003955,
                0.003468,
                0.003032,
                0.002597,
                0.001999,
                0.001401,
                0.000952,
                0.000503,
                0.000408,
                0.000314,
                0.000278,
                0.000242,
                0.000278,
                0.000315,
                0.000296,
                0.000277,
                0.000258,
                0.000239,
                0.000221,
                0.000204,
                0.000169,
                0.000135,
                0.000101,
                0.000067,
                0.000151,
                0.000236,
                0.000219,
                0.000202,
                0.000271,
                0.000340,
                0.000219,
                0.000099,
                0.000133,
                0.000167,
                0.000132,
                0.000098,
                0.000098,
                0.000098,
                0.000131,
                0.000165,
                0.000148,
                0.000131,
                0.000097,
                0.000064,
                0.000064,
                0.000064,
                0.000115,
                0.000165,
                0.000197,
                0.000229,
                0.000196,
                0.000164,
                0.000097,
                0.000031,
                0.000015,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000015,
                0.000030,
                0.000015,
                0.000000,
                0.000000,
                0.000000,
                0.000065,
                0.000131,
                0.000065,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000136,
                0.000271,
                0.000136,
                0.000000,
                0.000000,
                0.000000,
                0.000053,
                0.000106,
                0.000107,
                0.000107,
                0.000144,
                0.000180,
                0.000126,
                0.000071,
                0.000164,
                0.000258,
                0.000184,
                0.000111,
                0.000148,
                0.000186,
                0.000185,
                0.000185,
                0.000128,
                0.000071,
                0.000053,
                0.000035,
                0.000094,
                0.000153,
                0.000218,
                0.000282,
                0.000181,
                0.000079,
                0.000162,
                0.000245,
                0.000123,
                0.000000,
                0.000076,
                0.000151,
                0.000147,
                0.000143,
                0.000072,
                0.000000,
                0.000116,
                0.000231,
                0.000211,
                0.000190,
                0.000140,
                0.000091,
                0.000178,
                0.000266,
                0.000218,
                0.000170,
                0.000140,
                0.000111,
                0.000208,
                0.000304,
                0.000152,
                0.000000,
                0.000026,
                0.000052,
                0.000053,
                0.000054,
                0.000027,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000057,
                0.000114,
                0.000084,
                0.000055,
                0.000027,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000085,
                0.000170,
                0.000085,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000053,
                0.000107,
                0.000053,
                0.000000,
                0.000000,
                0.000000,
                0.000034,
                0.000069,
                0.000051,
                0.000034,
                0.000017,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000082,
                0.000164,
                0.000082,
                0.000000,
                0.000020,
                0.000040,
                0.000020,
                0.000000,
                0.000000,
                0.000000,
                0.000087,
                0.000174,
                0.000087,
                0.000000,
                0.000021,
                0.000042,
                0.000021,
                0.000000,
                0.000114,
                0.000227,
                0.000114,
                0.000000,
                0.000140,
                0.000280,
                0.000140,
                0.000000,
                0.000071,
                0.000142,
                0.000169,
                0.000195,
                0.000145,
                0.000095,
                0.000047,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000051,
                0.000101,
                0.000101,
                0.000101,
                0.000154,
                0.000207,
                0.000257,
                0.000306,
                0.000203,
                0.000100,
                0.000153,
                0.000206,
                0.000236,
                0.000266,
                0.000133,
                0.000000,
                0.000237,
                0.000475,
                0.000262,
                0.000050,
                0.000132,
                0.000215,
                0.000214,
                0.000213,
                0.000183,
                0.000153,
                0.000128,
                0.000103,
                0.000076,
                0.000050,
                0.000025,
                0.000000,
                0.000108,
                0.000216,
                0.000161,
                0.000105,
                0.000053,
                0.000000,
                0.000105,
                0.000210,
                0.000293,
                0.000377,
                0.000268,
                0.000159,
                0.000326,
                0.000493,
                0.000408,
                0.000323,
                0.000561,
                0.000798,
                0.000693,
                0.000587,
                0.000474,
                0.000361,
                0.000631,
                0.000900,
                0.000800,
                0.000699,
                0.000401,
                0.000102,
                0.000157,
                0.000212,
                0.000182,
                0.000152,
                0.000283,
                0.000414,
                0.000366,
                0.000319,
                0.000369,
                0.000419,
                0.000343,
                0.000267,
                0.000299,
                0.000330,
                0.000165,
                0.000000,
                0.000246,
                0.000492,
                0.000615,
                0.000738,
                0.000933,
                0.001128,
                0.000936,
                0.000744,
                0.000825,
                0.000906,
                0.001069,
                0.001233,
                0.001487,
                0.001741,
                0.001772,
                0.001803,
                0.001685,
                0.001567,
                0.001628,
                0.001689,
                0.001878,
                0.002067,
                0.002089,
                0.002112,
                0.002273,
                0.002434,
                0.002660,
                0.002886,
                0.003077,
                0.003269,
                0.003276,
                0.003283,
                0.003490,
                0.003697,
                0.003748,
                0.003799,
                0.003796,
                0.003793,
                0.003751,
                0.003710,
                0.003711,
                0.003712,
                0.003723,
                0.003734,
                0.004066,
                0.004397,
                0.004223,
                0.004049,
                0.004704,
                0.005358,
                0.005138,
                0.004919,
                0.004936,
                0.004954,
                0.005029,
                0.005104,
                0.005118,
                0.005133,
                0.005155,
                0.005178,
                0.005278,
                0.005378,
                0.005460,
                0.005542,
                0.005556,
                0.005570,
                0.005453,
                0.005336,
                0.005502,
                0.005668,
                0.005778,
                0.005887,
                0.005931,
                0.005975,
                0.006297,
                0.006619,
                0.006667,
                0.006715,
                0.006949,
                0.007184,
                0.007349,
                0.007514,
                0.008080,
                0.008647,
                0.008924,
                0.009201,
                0.009568,
                0.009934,
                0.010568,
                0.011202,
                0.011724,
                0.012246,
                0.012158,
                0.012070,
                0.006035,
                0.000000,
            ]
            self.specBandInfo["BAND5"] = specBandObj

            # BAND6
            specBandObj.phyBandName = "BAND6"
            specBandObj.bandID = "1"
            specBandObj.imgRes = 20
            specBandObj.wvLenMin = 360  # 520
            specBandObj.wvLenMax = 1000  # 590
            specBandObj.wvLenCen = 555.0
            specBandObj.respFuncStep = 1
            # Values imported from TerraLib
            specBandObj.respFunc = [
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000209,
                0.000418,
                0.000209,
                0.000000,
                0.000195,
                0.000389,
                0.000267,
                0.000144,
                0.000573,
                0.001002,
                0.000671,
                0.000339,
                0.000169,
                0.000000,
                0.000250,
                0.000500,
                0.000250,
                0.000000,
                0.000584,
                0.001167,
                0.000584,
                0.000000,
                0.000293,
                0.000586,
                0.000342,
                0.000098,
                0.000170,
                0.000241,
                0.000526,
                0.000811,
                0.000449,
                0.000087,
                0.000612,
                0.001138,
                0.000867,
                0.000596,
                0.001023,
                0.001450,
                0.001677,
                0.001904,
                0.001879,
                0.001853,
                0.002122,
                0.002390,
                0.002641,
                0.002892,
                0.002632,
                0.002372,
                0.002403,
                0.002435,
                0.002235,
                0.002035,
                0.001556,
                0.001078,
                0.001159,
                0.001240,
                0.001411,
                0.001582,
                0.001377,
                0.001173,
                0.000929,
                0.000685,
                0.000895,
                0.001105,
                0.000552,
                0.000000,
                0.000476,
                0.000951,
                0.001019,
                0.001086,
                0.000801,
                0.000515,
                0.001008,
                0.001501,
                0.001414,
                0.001327,
                0.000762,
                0.000196,
                0.000445,
                0.000694,
                0.000970,
                0.001247,
                0.001341,
                0.001435,
                0.001762,
                0.002089,
                0.002505,
                0.002921,
                0.002901,
                0.002881,
                0.003862,
                0.004842,
                0.005480,
                0.006117,
                0.008580,
                0.011043,
                0.010900,
                0.010757,
                0.010204,
                0.009650,
                0.010261,
                0.010873,
                0.014424,
                0.017975,
                0.028313,
                0.038651,
                0.058785,
                0.078918,
                0.107101,
                0.135284,
                0.175280,
                0.215276,
                0.282255,
                0.349233,
                0.443724,
                0.538215,
                0.617065,
                0.695915,
                0.728827,
                0.761740,
                0.759479,
                0.757218,
                0.761411,
                0.765605,
                0.765936,
                0.766268,
                0.779257,
                0.792247,
                0.805359,
                0.818471,
                0.827680,
                0.836890,
                0.831630,
                0.826370,
                0.830763,
                0.835156,
                0.838075,
                0.840994,
                0.851884,
                0.862775,
                0.871613,
                0.880452,
                0.880241,
                0.880030,
                0.885484,
                0.890939,
                0.885211,
                0.879483,
                0.877130,
                0.874777,
                0.886654,
                0.898532,
                0.908330,
                0.918129,
                0.929224,
                0.940320,
                0.943460,
                0.946601,
                0.951758,
                0.956916,
                0.962154,
                0.967392,
                0.965531,
                0.963671,
                0.965432,
                0.967193,
                0.972763,
                0.978333,
                0.977361,
                0.976390,
                0.976102,
                0.975815,
                0.976892,
                0.977968,
                0.983664,
                0.989360,
                0.993406,
                0.997452,
                0.998726,
                1.0,
                0.998557,
                0.997113,
                0.992439,
                0.987765,
                0.972028,
                0.956290,
                0.931981,
                0.907672,
                0.881162,
                0.854653,
                0.807615,
                0.760577,
                0.690585,
                0.620593,
                0.537942,
                0.455291,
                0.389490,
                0.323690,
                0.267297,
                0.210905,
                0.172901,
                0.134897,
                0.106833,
                0.078768,
                0.062850,
                0.046933,
                0.037085,
                0.027236,
                0.022045,
                0.016853,
                0.014417,
                0.011981,
                0.009940,
                0.007899,
                0.006773,
                0.005647,
                0.004798,
                0.003950,
                0.003830,
                0.003710,
                0.003684,
                0.003659,
                0.003680,
                0.003701,
                0.003553,
                0.003406,
                0.003418,
                0.003431,
                0.003290,
                0.003148,
                0.002717,
                0.002287,
                0.002099,
                0.001910,
                0.001597,
                0.001284,
                0.001096,
                0.000907,
                0.000763,
                0.000619,
                0.000600,
                0.000582,
                0.000728,
                0.000875,
                0.000726,
                0.000576,
                0.000645,
                0.000714,
                0.000548,
                0.000382,
                0.000516,
                0.000649,
                0.000610,
                0.000570,
                0.000285,
                0.000000,
                0.000126,
                0.000252,
                0.000410,
                0.000568,
                0.000630,
                0.000692,
                0.000770,
                0.000847,
                0.000532,
                0.000216,
                0.000276,
                0.000335,
                0.000329,
                0.000322,
                0.000276,
                0.000231,
                0.000363,
                0.000495,
                0.000335,
                0.000175,
                0.000304,
                0.000433,
                0.000375,
                0.000316,
                0.000243,
                0.000170,
                0.000085,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000019,
                0.000037,
                0.000128,
                0.000218,
                0.000190,
                0.000161,
                0.000081,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000113,
                0.000226,
                0.000166,
                0.000105,
                0.000362,
                0.000619,
                0.000723,
                0.000827,
                0.000858,
                0.000889,
                0.002263,
                0.003638,
                0.004899,
                0.006161,
                0.006254,
                0.006347,
                0.005262,
                0.004178,
                0.002662,
                0.001146,
                0.000573,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000059,
                0.000117,
                0.000059,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000125,
                0.000251,
                0.000218,
                0.000186,
                0.000093,
                0.000000,
                0.000266,
                0.000532,
                0.000434,
                0.000336,
                0.000914,
                0.001491,
                0.001188,
                0.000885,
                0.000747,
                0.000609,
                0.000363,
                0.000117,
                0.000059,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000056,
                0.000112,
                0.000056,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000023,
                0.000045,
                0.000023,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
            ]
            self.specBandInfo["BAND6"] = specBandObj

            # BAND7
            specBandObj.phyBandName = "BAND7"
            specBandObj.bandID = "2"
            specBandObj.imgRes = 20
            specBandObj.wvLenMin = 360  # 630
            specBandObj.wvLenMax = 1000  # 690
            specBandObj.wvLenCen = 660.0
            specBandObj.respFuncStep = 1
            # Values imported from TerraLib
            specBandObj.respFunc = [
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000222,
                0.000445,
                0.000222,
                0.000000,
                0.000320,
                0.000639,
                0.000320,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000380,
                0.000761,
                0.000380,
                0.000000,
                0.000146,
                0.000291,
                0.000146,
                0.000000,
                0.000042,
                0.000084,
                0.000262,
                0.000440,
                0.000220,
                0.000000,
                0.001105,
                0.002210,
                0.001105,
                0.000000,
                0.000000,
                0.000000,
                0.000033,
                0.000067,
                0.000033,
                0.000000,
                0.000165,
                0.000329,
                0.000165,
                0.000000,
                0.000219,
                0.000439,
                0.000372,
                0.000305,
                0.000240,
                0.000176,
                0.000292,
                0.000408,
                0.000204,
                0.000000,
                0.000360,
                0.000721,
                0.000495,
                0.000269,
                0.000135,
                0.000000,
                0.000282,
                0.000563,
                0.000282,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.001060,
                0.002120,
                0.002134,
                0.002148,
                0.001074,
                0.000000,
                0.000000,
                0.000000,
                0.000247,
                0.000495,
                0.000378,
                0.000261,
                0.000131,
                0.000000,
                0.000124,
                0.000249,
                0.000211,
                0.000174,
                0.000311,
                0.000448,
                0.000342,
                0.000236,
                0.000268,
                0.000300,
                0.000462,
                0.000625,
                0.000425,
                0.000225,
                0.000380,
                0.000536,
                0.000502,
                0.000468,
                0.000371,
                0.000274,
                0.000453,
                0.000631,
                0.000506,
                0.000382,
                0.000585,
                0.000788,
                0.000666,
                0.000544,
                0.000313,
                0.000082,
                0.000307,
                0.000532,
                0.000306,
                0.000081,
                0.000323,
                0.000566,
                0.000295,
                0.000025,
                0.000135,
                0.000245,
                0.000296,
                0.000348,
                0.000425,
                0.000503,
                0.000678,
                0.000853,
                0.000822,
                0.000791,
                0.000664,
                0.000537,
                0.000611,
                0.000686,
                0.000954,
                0.001222,
                0.001048,
                0.000874,
                0.000768,
                0.000661,
                0.000811,
                0.000961,
                0.000908,
                0.000855,
                0.000655,
                0.000455,
                0.000311,
                0.000167,
                0.000542,
                0.000916,
                0.000796,
                0.000676,
                0.000514,
                0.000353,
                0.000445,
                0.000537,
                0.000649,
                0.000761,
                0.000737,
                0.000712,
                0.000781,
                0.000850,
                0.000738,
                0.000627,
                0.000553,
                0.000478,
                0.000409,
                0.000340,
                0.000524,
                0.000708,
                0.000843,
                0.000979,
                0.000909,
                0.000840,
                0.000812,
                0.000785,
                0.000740,
                0.000695,
                0.000672,
                0.000649,
                0.000849,
                0.001050,
                0.000984,
                0.000917,
                0.000986,
                0.001054,
                0.000852,
                0.000649,
                0.000689,
                0.000729,
                0.000731,
                0.000734,
                0.000760,
                0.000786,
                0.000983,
                0.001181,
                0.001199,
                0.001217,
                0.001443,
                0.001669,
                0.002011,
                0.002352,
                0.002776,
                0.003199,
                0.003541,
                0.003882,
                0.005590,
                0.007298,
                0.008868,
                0.010439,
                0.012869,
                0.015298,
                0.019036,
                0.022773,
                0.029037,
                0.035301,
                0.048376,
                0.061450,
                0.085476,
                0.109502,
                0.157513,
                0.205525,
                0.282238,
                0.358951,
                0.454161,
                0.549372,
                0.624941,
                0.700511,
                0.751429,
                0.802346,
                0.837077,
                0.871807,
                0.890203,
                0.908599,
                0.918331,
                0.928063,
                0.928965,
                0.929867,
                0.934336,
                0.938805,
                0.939611,
                0.940417,
                0.949866,
                0.959316,
                0.955100,
                0.950884,
                0.949815,
                0.948746,
                0.941605,
                0.934463,
                0.948301,
                0.962139,
                0.972099,
                0.982059,
                0.991029,
                1.0,
                0.988619,
                0.977237,
                0.973996,
                0.970755,
                0.968356,
                0.965958,
                0.955100,
                0.944241,
                0.935509,
                0.926778,
                0.927118,
                0.927459,
                0.918239,
                0.909019,
                0.908467,
                0.907914,
                0.907296,
                0.906678,
                0.900698,
                0.894718,
                0.881502,
                0.868286,
                0.859020,
                0.849755,
                0.811625,
                0.773495,
                0.694240,
                0.614984,
                0.498664,
                0.382344,
                0.297868,
                0.213392,
                0.157202,
                0.101012,
                0.077027,
                0.053043,
                0.041347,
                0.029652,
                0.023894,
                0.018137,
                0.015042,
                0.011948,
                0.010298,
                0.008648,
                0.007603,
                0.006558,
                0.006376,
                0.006193,
                0.005469,
                0.004745,
                0.004184,
                0.003622,
                0.003565,
                0.003508,
                0.003398,
                0.003289,
                0.003129,
                0.002969,
                0.002914,
                0.002858,
                0.002438,
                0.002017,
                0.001589,
                0.001160,
                0.001230,
                0.001300,
                0.001207,
                0.001114,
                0.000926,
                0.000737,
                0.000833,
                0.000929,
                0.000819,
                0.000709,
                0.000759,
                0.000810,
                0.000810,
                0.000811,
                0.000835,
                0.000860,
                0.000922,
                0.000985,
                0.001021,
                0.001056,
                0.001163,
                0.001270,
                0.001196,
                0.001123,
                0.001090,
                0.001057,
                0.001137,
                0.001217,
                0.001141,
                0.001064,
                0.001231,
                0.001398,
                0.001169,
                0.000941,
                0.001224,
                0.001507,
                0.001102,
                0.000697,
                0.000908,
                0.001118,
                0.001046,
                0.000973,
                0.000783,
                0.000592,
                0.000413,
                0.000233,
                0.000616,
                0.000999,
                0.001172,
                0.001345,
                0.001323,
                0.001302,
                0.001324,
                0.001346,
                0.001235,
                0.001124,
                0.001178,
                0.001233,
                0.001055,
                0.000878,
                0.001235,
                0.001592,
                0.001512,
                0.001432,
                0.001537,
                0.001641,
                0.001446,
                0.001252,
                0.001105,
                0.000959,
                0.001041,
                0.001122,
                0.001257,
                0.001392,
                0.001199,
                0.001006,
                0.001140,
                0.001274,
                0.001402,
                0.001529,
                0.001460,
                0.001391,
                0.001443,
                0.001495,
                0.001428,
                0.001361,
                0.001640,
                0.001918,
                0.001792,
                0.001667,
                0.001649,
                0.001631,
                0.001429,
                0.001226,
                0.000928,
                0.000630,
                0.000570,
                0.000511,
                0.000444,
                0.000377,
                0.000929,
                0.001482,
                0.001459,
                0.001435,
                0.001143,
                0.000851,
                0.000954,
                0.001058,
                0.001149,
                0.001240,
                0.000813,
                0.000385,
                0.000563,
                0.000741,
                0.000906,
                0.001070,
                0.000831,
                0.000592,
                0.000488,
                0.000385,
                0.000426,
                0.000467,
                0.000234,
                0.000000,
                0.000764,
                0.001528,
                0.001219,
                0.000909,
                0.000948,
                0.000988,
                0.001165,
                0.001342,
                0.000757,
                0.000172,
                0.000569,
                0.000966,
                0.000902,
                0.000837,
                0.000418,
                0.000000,
                0.000793,
                0.001586,
                0.001178,
                0.000770,
                0.000615,
                0.000461,
                0.000570,
                0.000680,
                0.000466,
                0.000252,
                0.000434,
                0.000616,
                0.000436,
                0.000256,
                0.000473,
                0.000691,
                0.000899,
                0.001106,
                0.001038,
                0.000969,
                0.000954,
                0.000940,
                0.000771,
                0.000602,
                0.001078,
                0.001554,
                0.001260,
                0.000965,
                0.001007,
                0.001049,
                0.000751,
                0.000452,
                0.000942,
                0.001432,
                0.001489,
                0.001547,
                0.001356,
                0.001165,
                0.001582,
                0.001999,
                0.001913,
                0.001827,
                0.001480,
                0.001133,
                0.001432,
                0.001731,
                0.001272,
                0.000813,
                0.001042,
                0.001271,
                0.001188,
                0.001105,
                0.001142,
                0.001179,
                0.000678,
                0.000177,
                0.000543,
                0.000909,
                0.000911,
                0.000913,
                0.001125,
                0.001338,
                0.001307,
                0.001275,
                0.001404,
                0.001532,
                0.001218,
                0.000904,
                0.001173,
                0.001442,
                0.001307,
                0.001172,
                0.001640,
                0.002109,
                0.001465,
                0.000821,
                0.001337,
                0.001853,
                0.001729,
                0.001605,
                0.001342,
                0.001080,
                0.000776,
                0.000472,
                0.000438,
                0.000405,
                0.000550,
                0.000695,
                0.000769,
                0.000843,
                0.000914,
                0.000986,
                0.001044,
                0.001103,
                0.001157,
                0.001211,
                0.001170,
                0.001128,
                0.000853,
                0.000579,
                0.000289,
                0.000000,
                0.000602,
                0.001204,
                0.001200,
                0.001196,
                0.001495,
                0.001794,
                0.001541,
                0.001288,
                0.000817,
                0.000346,
                0.000582,
                0.000818,
                0.000503,
                0.000188,
                0.000265,
                0.000343,
                0.000535,
                0.000728,
                0.000842,
                0.000957,
                0.000692,
                0.000426,
                0.000431,
                0.000435,
                0.000278,
                0.000121,
                0.000284,
                0.000448,
                0.000362,
                0.000277,
                0.000423,
                0.000569,
                0.000615,
                0.000660,
                0.000553,
                0.000445,
                0.000222,
                0.000000,
            ]
            self.specBandInfo["BAND7"] = specBandObj

            # BAND8
            specBandObj.phyBandName = "BAND8"
            specBandObj.bandID = "3"
            specBandObj.imgRes = 20
            specBandObj.wvLenMin = 360  # 770
            specBandObj.wvLenMax = 1000  # 890
            specBandObj.wvLenCen = 830.0
            specBandObj.respFuncStep = 1
            # Values imported from TerraLib
            specBandObj.respFunc = [
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000090,
                0.000181,
                0.000090,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000064,
                0.000129,
                0.000064,
                0.000000,
                0.000126,
                0.000253,
                0.000126,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000103,
                0.000205,
                0.000103,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000033,
                0.000066,
                0.000033,
                0.000000,
                0.000000,
                0.000000,
                0.000063,
                0.000127,
                0.000063,
                0.000000,
                0.000000,
                0.000000,
                0.000027,
                0.000055,
                0.000027,
                0.000000,
                0.000000,
                0.000000,
                0.000056,
                0.000112,
                0.000081,
                0.000051,
                0.000025,
                0.000000,
                0.000052,
                0.000105,
                0.000185,
                0.000266,
                0.000556,
                0.000847,
                0.000447,
                0.000047,
                0.000046,
                0.000046,
                0.000046,
                0.000045,
                0.000023,
                0.000000,
                0.000046,
                0.000093,
                0.000068,
                0.000043,
                0.000021,
                0.000000,
                0.000000,
                0.000000,
                0.000020,
                0.000041,
                0.000020,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000018,
                0.000037,
                0.000018,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000017,
                0.000034,
                0.000017,
                0.000000,
                0.000017,
                0.000033,
                0.000017,
                0.000000,
                0.000016,
                0.000033,
                0.000016,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000016,
                0.000032,
                0.000032,
                0.000032,
                0.000016,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000016,
                0.000031,
                0.000016,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000015,
                0.000031,
                0.000031,
                0.000031,
                0.000015,
                0.000000,
                0.000000,
                0.000000,
                0.000033,
                0.000065,
                0.000033,
                0.000000,
                0.000000,
                0.000000,
                0.000015,
                0.000031,
                0.000015,
                0.000000,
                0.000000,
                0.000000,
                0.000016,
                0.000031,
                0.000031,
                0.000031,
                0.000015,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000015,
                0.000031,
                0.000048,
                0.000066,
                0.000033,
                0.000000,
                0.000016,
                0.000031,
                0.000032,
                0.000032,
                0.000069,
                0.000106,
                0.000145,
                0.000184,
                0.000166,
                0.000149,
                0.000149,
                0.000149,
                0.000149,
                0.000149,
                0.000150,
                0.000151,
                0.000171,
                0.000191,
                0.000153,
                0.000115,
                0.000115,
                0.000114,
                0.000114,
                0.000114,
                0.000094,
                0.000073,
                0.000134,
                0.000195,
                0.000198,
                0.000200,
                0.000162,
                0.000123,
                0.000276,
                0.000429,
                0.000364,
                0.000299,
                0.000211,
                0.000123,
                0.000282,
                0.000440,
                0.000371,
                0.000302,
                0.000204,
                0.000107,
                0.000139,
                0.000171,
                0.000201,
                0.000231,
                0.000261,
                0.000290,
                0.000221,
                0.000153,
                0.000165,
                0.000177,
                0.000235,
                0.000293,
                0.000320,
                0.000347,
                0.000245,
                0.000143,
                0.000112,
                0.000082,
                0.000054,
                0.000026,
                0.000013,
                0.000000,
                0.000000,
                0.000000,
                0.000013,
                0.000025,
                0.000026,
                0.000026,
                0.000013,
                0.000000,
                0.000000,
                0.000000,
                0.000014,
                0.000027,
                0.000027,
                0.000027,
                0.000058,
                0.000088,
                0.000072,
                0.000057,
                0.000057,
                0.000058,
                0.000043,
                0.000027,
                0.000058,
                0.000089,
                0.000074,
                0.000059,
                0.000044,
                0.000028,
                0.000076,
                0.000124,
                0.000092,
                0.000061,
                0.000062,
                0.000062,
                0.000095,
                0.000128,
                0.000147,
                0.000166,
                0.000300,
                0.000434,
                0.000832,
                0.001231,
                0.001661,
                0.002091,
                0.002649,
                0.003208,
                0.003317,
                0.003425,
                0.003841,
                0.004257,
                0.005741,
                0.007225,
                0.008940,
                0.010654,
                0.012552,
                0.014451,
                0.016514,
                0.018577,
                0.020399,
                0.022222,
                0.024722,
                0.027222,
                0.030597,
                0.033971,
                0.038271,
                0.042570,
                0.047825,
                0.053081,
                0.062215,
                0.071349,
                0.085534,
                0.099720,
                0.123384,
                0.147048,
                0.183884,
                0.220721,
                0.282534,
                0.344347,
                0.416241,
                0.488134,
                0.582517,
                0.676900,
                0.772149,
                0.867399,
                0.904229,
                0.941058,
                0.970529,
                1.0,
                0.968928,
                0.937855,
                0.925289,
                0.912723,
                0.875339,
                0.837955,
                0.836936,
                0.835917,
                0.828990,
                0.822063,
                0.814755,
                0.807447,
                0.799989,
                0.792531,
                0.797652,
                0.802773,
                0.813029,
                0.823285,
                0.826962,
                0.830639,
                0.828741,
                0.826842,
                0.818331,
                0.809821,
                0.812974,
                0.816128,
                0.810388,
                0.804647,
                0.811876,
                0.819104,
                0.813334,
                0.807563,
                0.805787,
                0.804011,
                0.796793,
                0.789574,
                0.784041,
                0.778508,
                0.781810,
                0.785112,
                0.782458,
                0.779803,
                0.770303,
                0.760804,
                0.761772,
                0.762740,
                0.746772,
                0.730804,
                0.730303,
                0.729802,
                0.736810,
                0.743817,
                0.738771,
                0.733725,
                0.729064,
                0.724403,
                0.711984,
                0.699564,
                0.700354,
                0.701145,
                0.701163,
                0.701182,
                0.705485,
                0.709787,
                0.711417,
                0.713047,
                0.698679,
                0.684310,
                0.690605,
                0.696899,
                0.695887,
                0.694874,
                0.689262,
                0.683651,
                0.669994,
                0.656338,
                0.660359,
                0.664380,
                0.667846,
                0.671312,
                0.655698,
                0.640085,
                0.646025,
                0.651965,
                0.645562,
                0.639159,
                0.624486,
                0.609814,
                0.602325,
                0.594837,
                0.592166,
                0.589496,
                0.581508,
                0.573521,
                0.574610,
                0.575700,
                0.567022,
                0.558344,
                0.546279,
                0.534213,
                0.526024,
                0.517834,
                0.492070,
                0.466307,
                0.447254,
                0.428202,
                0.392737,
                0.357272,
                0.323391,
                0.289510,
                0.247737,
                0.205963,
                0.177736,
                0.149509,
                0.126356,
                0.103203,
                0.088285,
                0.073367,
                0.064268,
                0.055169,
                0.045996,
                0.036822,
                0.032975,
                0.029129,
                0.025682,
                0.022235,
                0.018526,
                0.014817,
                0.013007,
                0.011198,
                0.009714,
                0.008230,
                0.006963,
                0.005696,
                0.005645,
                0.005593,
                0.005353,
                0.005112,
                0.004550,
                0.003987,
                0.003337,
                0.002686,
                0.001680,
                0.000674,
                0.000734,
                0.000793,
                0.000657,
                0.000520,
                0.000399,
                0.000277,
                0.000280,
                0.000283,
                0.000227,
                0.000171,
                0.000142,
                0.000113,
                0.000115,
                0.000116,
                0.000112,
                0.000107,
                0.000080,
                0.000052,
                0.000052,
                0.000052,
                0.000026,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000025,
                0.000051,
                0.000052,
                0.000053,
                0.000054,
                0.000054,
                0.000027,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000027,
                0.000053,
                0.000027,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
            ]
            self.specBandInfo["BAND8"] = specBandObj

            # Get 6S spectral response functions.
            for specRspBand in self.specBandInfo:
                respFuncSize = len(self.specBandInfo[specRspBand].respFunc)
                tmpWvLenMin = self.specBandInfo[specRspBand].wvLenMin
                tmpWvLenMax = self.specBandInfo[specRspBand].wvLenMax
                numWvLens = math.ceil(tmpWvLenMax - tmpWvLenMin)
                lenDiff = respFuncSize - numWvLens
                if lenDiff > 0:
                    tmpWvLenMax = tmpWvLenMax + lenDiff
                elif lenDiff < 0:
                    tmpWvLenMin = tmpWvLenMin + (lenDiff * (-1))
                wvLensIn = numpy.arange(tmpWvLenMin, tmpWvLenMax, 1)
                olWVLens, olSpecResp = arcsilib.arcsiutils.resampleSpectralResponseFunc(
                    wvLensIn, self.specBandInfo[specRspBand].respFunc, 2.5, "linear"
                )
                self.specBandInfo[specRspBand].respFunc6S = olSpecResp
                self.specBandInfo[specRspBand].wvLenMin6S = olWVLens[0] / 1000
                self.specBandInfo[specRspBand].wvLenMax6S = olWVLens[-1] / 1000

        except Exception as e:
            raise e


    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        raise ARCSIException("Not Implemented")

    # def getSolarIrrStdSolarGeom(self):
    #     """
    #     Get Solar Azimuth and Zenith as standard geometry.
    #     Azimuth: N=0, E=90, S=180, W=270.
    #     """
    #     solarAz = rsgislib.imagecalibration.solarangles.get_solar_irr_convention_solar_azimuth_from_usgs(
    #         self.solarAzimuth
    #     )
    #     return (solarAz, self.solarZenith)

    # def getSensorViewGeom(self):
    #     """
    #     Get sensor viewing angles
    #     returns (viewAzimuth, viewZenith)
    #     """
    #     return (0.0, 0.0)

    # def generateOutputBaseName(self):
    #     """
    #     Provides an implementation for the landsat sensor
    #     """
    #     rowpath = "r" + str(self.row) + "p" + str(self.path)
    #     outname = self.defaultGenBaseOutFileName()
    #     outname = outname + str("_") + rowpath
    #     return outname

    # def generateMetaDataFile(
    #     self,
    #     outputPath,
    #     outputFileName,
    #     productsStr,
    #     validMaskImage="",
    #     footprintCalc=False,
    #     calcdValuesDict=None,
    #     outFilesDict=None,
    # ):
    #     """
    #     Generate file metadata.
    #     """
    #     if outFilesDict is None:
    #         outFilesDict = dict()
    #     if calcdValuesDict is None:
    #         calcdValuesDict = dict()
    #     outJSONFilePath = os.path.join(outputPath, outputFileName)
    #     jsonData = self.getJSONDictDefaultMetaData(
    #         productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict
    #     )
    #     sensorInfo = jsonData["SensorInfo"]
    #     sensorInfo["Row"] = self.row
    #     sensorInfo["Path"] = self.path
    #     sensorInfo["SensorID"] = self.sensorID
    #     sensorInfo["SpacecraftID"] = self.spacecraftID
    #     acqDict = jsonData["AcquasitionInfo"]
    #     acqDict["EarthSunDistance"] = self.earthSunDistance
    #     imgInfo = dict()
    #     imgInfo["CloudCover"] = self.cloudCover
    #     imgInfo["CloudCoverLand"] = self.cloudCoverLand
    #     imgInfo["CellSizePan"] = self.gridCellSizePan
    #     imgInfo["CellSizeRefl"] = self.gridCellSizeRefl
    #     imgInfo["CellSizeTherm"] = self.gridCellSizeTherm
    #     jsonData["ImageInfo"] = imgInfo

    #     with open(outJSONFilePath, "w") as outfile:
    #         json.dump(
    #             jsonData,
    #             outfile,
    #             sort_keys=True,
    #             indent=4,
    #             separators=(",", ": "),
    #             ensure_ascii=False,
    #         )

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.band5File):
            imageDataPresent = False
        if not os.path.exists(self.band6File):
            imageDataPresent = False
        if not os.path.exists(self.band7File):
            imageDataPresent = False
        if not os.path.exists(self.band8File):
            imageDataPresent = False

        return imageDataPresent

    def hasThermal(self):
        return False

    def applyImageDataMask(
        self,
        inputHeader,
        inputImage,
        outputPath,
        outputMaskName,
        outputImgName,
        outFormat,
        outWKTFile,
    ):
        raise ARCSIException(
            "CBERS-4/MUX does not provide any image masks, do not use the MASK option."
        )

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(
        self,
        outputPath,
        resampleToLowResImg,
        resampleMethod=rsgislib.INTERP_CUBIC,
        multicore=False,
    ):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    # def generateValidImageDataMask(
    #     self, outputPath, outputMaskName, viewAngleImg, outFormat
    # ):
    #     print("Create the valid data mask")
    #     tmpBaseName = os.path.splitext(outputMaskName)[0]
    #     tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName + "vldpxlmsk.kea")
    #     outputImage = os.path.join(outputPath, outputMaskName)
    #     inImages = [
    #         self.band1File,
    #         self.band2File,
    #         self.band3File,
    #         self.band4File,
    #         self.band5File,
    #         self.band6File,
    #         self.band7File,
    #         self.band10File,
    #         self.band11File,
    #     ]
    #     rsgislib.imageutils.gen_valid_mask(
    #         input_imgs=inImages,
    #         output_img=tmpValidPxlMsk,
    #         gdalformat="KEA",
    #         no_data_val=0.0,
    #     )
    #     rsgislib.rastergis.pop_rat_img_stats(tmpValidPxlMsk, True, False, True)
    #     # Check there is valid data
    #     ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
    #     Histogram = rat.readColumn(ratDS, "Histogram")
    #     ratDS = None
    #     if Histogram.shape[0] < 2:
    #         raise ARCSIException("There is no valid data in this image.")
    #     if not os.path.exists(viewAngleImg):
    #         print("Calculate Image Angles.")
    #         imgInfo = rios.fileinfo.ImageInfo(tmpValidPxlMsk)
    #         corners = fmask.landsatangles.findImgCorners(tmpValidPxlMsk, imgInfo)
    #         nadirLine = fmask.landsatangles.findNadirLine(corners)
    #         extentSunAngles = fmask.landsatangles.sunAnglesForExtent(
    #             imgInfo, self.fmaskMTLInfo
    #         )
    #         satAzimuth = fmask.landsatangles.satAzLeftRight(nadirLine)
    #         fmask.landsatangles.makeAnglesImage(
    #             tmpValidPxlMsk,
    #             viewAngleImg,
    #             nadirLine,
    #             extentSunAngles,
    #             satAzimuth,
    #             imgInfo,
    #         )
    #         dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
    #         if not dataset is None:
    #             dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
    #             dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
    #             dataset.GetRasterBand(3).SetDescription("SolorAzimuth")
    #             dataset.GetRasterBand(4).SetDescription("SolorZenith")
    #         dataset = None
    #     rsgislib.imagecalc.band_math(
    #         outputImage,
    #         "(VA<14)&&(VM==1)?1:0",
    #         outFormat,
    #         rsgislib.TYPE_8UINT,
    #         [
    #             rsgislib.imagecalc.BandDefn("VA", viewAngleImg, 2),
    #             rsgislib.imagecalc.BandDefn("VM", tmpValidPxlMsk, 1),
    #         ],
    #     )

    #     rsgislib.imageutils.delete_gdal_layer(tmpValidPxlMsk)
    #     return outputImage

    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple(
            "LSBand", ["band_name", "input_img", "img_band", "add_val", "multi_val"]
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Blue",
                input_img=self.band5File,
                img_band=1,
                add_val=self.b5RadAdd,
                multi_val=self.b5RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band6File,
                img_band=1,
                add_val=self.b6RadAdd,
                multi_val=self.b6RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band7File,
                img_band=1,
                add_val=self.b7RadAdd,
                multi_val=self.b7RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR",
                input_img=self.band8File,
                img_band=1,
                add_val=self.b8RadAdd,
                multi_val=self.b8RadMulti,
            )
        )
        processed_bands = []
        profile = None
        for i, band_cfg in enumerate(bandDefnSeq):
            with rasterio.open(band_cfg['file']) as src:
                data = src.read(1).astype(numpy.float32)

                radiance_data = band_cfg['mult'] * data + band_cfg['add']

                processed_bands.append(radiance_data)

                # Guardar perfil da primeira banda
                if i == 0:
                    profile = src.profile.copy()
                    profile.update({
                        'count': len(bandDefnSeq),
                        'dtype': rasterio.float32,
                        'driver': outFormat,
                        'nodata': None
                    })
        
        # Empilhar todas as bandas
        if len(processed_bands) > 0:
            stacked_data = numpy.stack(processed_bands, axis=0)
            
            # Escrever arquivo de saída
            with rasterio.open(outputImage, 'w', **profile) as dst:
                dst.write(stacked_data)
        else:
            raise ValueError("Nenhuma banda processada")

        # rsgislib.imagecalibration.landsat_to_radiance_multi_add(
        #     outputImage, outFormat, bandDefnSeq
        # )

        return outputImage, None

    # def generateImageSaturationMask(self, outputPath, outputName, outFormat):
    #     print("Generate Saturation Image")
    #     outputImage = os.path.join(outputPath, outputName)

    #     lsBand = collections.namedtuple(
    #         "LSBand", ["band_name", "input_img", "img_band", "sat_val"]
    #     )
    #     bandDefnSeq = list()
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="Coastal",
    #             input_img=self.band1File,
    #             img_band=1,
    #             sat_val=self.b1CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="Blue",
    #             input_img=self.band2File,
    #             img_band=1,
    #             sat_val=self.b2CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="Green",
    #             input_img=self.band3File,
    #             img_band=1,
    #             sat_val=self.b3CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="Red",
    #             input_img=self.band4File,
    #             img_band=1,
    #             sat_val=self.b4CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="NIR",
    #             input_img=self.band5File,
    #             img_band=1,
    #             sat_val=self.b5CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="SWIR1",
    #             input_img=self.band6File,
    #             img_band=1,
    #             sat_val=self.b6CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="SWIR2",
    #             input_img=self.band7File,
    #             img_band=1,
    #             sat_val=self.b7CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="ThermalB10",
    #             input_img=self.band10File,
    #             img_band=1,
    #             sat_val=self.b10CalMax,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="ThermalB11",
    #             input_img=self.band11File,
    #             img_band=1,
    #             sat_val=self.b11CalMax,
    #         )
    #     )

    #     rsgislib.imagecalibration.saturated_pixels_mask(
    #         outputImage, outFormat, bandDefnSeq
    #     )

    #     return outputImage

    #     lsBand = collections.namedtuple("LSBand", ["band_name", "img_band", "k1", "k2"])
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="ThermalB10",
    #             img_band=1,
    #             k1=self.k1ConstB10,
    #             k2=self.k2ConstB10,
    #         )
    #     )
    #     bandDefnSeq.append(
    #         lsBand(
    #             band_name="ThermalB11",
    #             img_band=2,
    #             k1=self.k1ConstB11,
    #             k2=self.k2ConstB11,
    #         )
    #     )
    #     rsgislib.imagecalibration.landsat_thermal_rad_to_brightness(
    #         inputRadImage,
    #         outputThermalImage,
    #         outFormat,
    #         rsgislib.TYPE_32INT,
    #         scaleFactor,
    #         bandDefnSeq,
    #     )
    #     return outputThermalImage

    # def convertImageToTOARefl(
    #     self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    # ):
    #     print("Converting to TOA")
    #     outputImage = os.path.join(outputPath, outputName)
    #     solarIrradianceVals = list()
    #     IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
    #     solarIrradianceVals.append(IrrVal(irradiance=1876.61))
    #     solarIrradianceVals.append(IrrVal(irradiance=1970.03))
    #     solarIrradianceVals.append(IrrVal(irradiance=1848.9))
    #     solarIrradianceVals.append(IrrVal(irradiance=1571.3))
    #     solarIrradianceVals.append(IrrVal(irradiance=967.66))
    #     solarIrradianceVals.append(IrrVal(irradiance=245.73))
    #     solarIrradianceVals.append(IrrVal(irradiance=82.03))
    #     rsgislib.imagecalibration.radiance_to_toa_refl(
    #         inputRadImage,
    #         outputImage,
    #         outFormat,
    #         rsgislib.TYPE_16UINT,
    #         scaleFactor,
    #         self.acquisitionTime.year,
    #         self.acquisitionTime.month,
    #         self.acquisitionTime.day,
    #         self.solarZenith,
    #         solarIrradianceVals,
    #     )
    #     return outputImage

    # def generateCloudMask(
    #     self,
    #     inputReflImage,
    #     inputSatImage,
    #     inputThermalImage,
    #     inputViewAngleImg,
    #     inputValidImg,
    #     outputPath,
    #     outputCloudName,
    #     outputCloudProb,
    #     outFormat,
    #     tmpPath,
    #     scaleFactor,
    #     cloud_msk_methods=None,
    # ):
    #     import rsgislib.imageutils

    #     try:
    #         outputImage = os.path.join(outputPath, outputCloudName)
    #         tmpBaseName = os.path.splitext(outputCloudName)[0]
    #         tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

    #         tmpDIRExisted = True
    #         if not os.path.exists(tmpBaseDIR):
    #             os.makedirs(tmpBaseDIR)
    #             tmpDIRExisted = False

    #         if (cloud_msk_methods is None) or (cloud_msk_methods == "FMASK"):
    #             tmpFMaskOut = os.path.join(tmpBaseDIR, tmpBaseName + "_pyfmaskout.kea")

    #             # Create tmp TOA stack with Band 9 (Cirrus)
    #             tmpTOAB9Out = os.path.join(tmpBaseDIR, tmpBaseName + "_TOAB9.kea")
    #             toaExp = (
    #                 "((b1*"
    #                 + str(self.b9ReflMulti)
    #                 + ")+"
    #                 + str(self.b9ReflAdd)
    #                 + ")*"
    #                 + str(scaleFactor)
    #             )
    #             rsgislib.imagecalc.image_math(
    #                 self.band9File,
    #                 tmpTOAB9Out,
    #                 toaExp,
    #                 "KEA",
    #                 rsgislib.TYPE_16UINT,
    #                 False,
    #             )
    #             if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
    #                 inputReflImage, tmpTOAB9Out
    #             ):
    #                 tmpTOAB9OutNotProj = tmpTOAB9Out
    #                 tmpTOAB9Out = os.path.join(
    #                     tmpBaseDIR, tmpBaseName + "_TOAB9_reproj.kea"
    #                 )
    #                 rsgislib.imageutils.resample_img_to_match(
    #                     inputReflImage,
    #                     tmpTOAB9OutNotProj,
    #                     tmpTOAB9Out,
    #                     "KEA",
    #                     rsgislib.INTERP_CUBIC,
    #                     rsgislib.TYPE_16UINT,
    #                 )

    #             tmpReflStackOut = os.path.join(
    #                 tmpBaseDIR, tmpBaseName + "_TOAreflStack.kea"
    #             )
    #             rsgislib.imageutils.stack_img_bands(
    #                 [inputReflImage, tmpTOAB9Out],
    #                 None,
    #                 tmpReflStackOut,
    #                 None,
    #                 0,
    #                 "KEA",
    #                 rsgislib.TYPE_16UINT,
    #             )

    #             tmpThermStack = os.path.join(
    #                 tmpBaseDIR, tmpBaseName + "_ThermB10B11Stack.kea"
    #             )
    #             rsgislib.imageutils.stack_img_bands(
    #                 [self.band10File, self.band11File],
    #                 None,
    #                 tmpThermStack,
    #                 None,
    #                 0,
    #                 "KEA",
    #                 rsgislib.TYPE_16UINT,
    #             )

    #             tmpThermalLayer = tmpThermStack
    #             if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
    #                 inputThermalImage, tmpThermStack
    #             ):
    #                 tmpThermalLayer = os.path.join(
    #                     tmpBaseDIR, tmpBaseName + "_thermalresample.kea"
    #                 )
    #                 rsgislib.imageutils.resample_img_to_match(
    #                     inputThermalImage,
    #                     tmpThermStack,
    #                     tmpThermalLayer,
    #                     "KEA",
    #                     rsgislib.INTERP_CUBIC,
    #                     rsgislib.TYPE_32FLOAT,
    #                 )

    #             minCloudSize = 0
    #             cloudBufferDistance = 150
    #             shadowBufferDistance = 300

    #             fmaskFilenames = fmask.config.FmaskFilenames()
    #             fmaskFilenames.setTOAReflectanceFile(tmpReflStackOut)
    #             fmaskFilenames.setThermalFile(tmpThermalLayer)
    #             fmaskFilenames.setSaturationMask(inputSatImage)
    #             fmaskFilenames.setOutputCloudMaskFile(tmpFMaskOut)

    #             thermalGain1040um = self.b10RadMulti
    #             thermalOffset1040um = self.b10RadAdd
    #             thermalBand1040um = 0
    #             thermalInfo = fmask.config.ThermalFileInfo(
    #                 thermalBand1040um,
    #                 thermalGain1040um,
    #                 thermalOffset1040um,
    #                 self.k1ConstB10,
    #                 self.k2ConstB10,
    #             )

    #             anglesInfo = fmask.config.AnglesFileInfo(
    #                 inputViewAngleImg,
    #                 3,
    #                 inputViewAngleImg,
    #                 2,
    #                 inputViewAngleImg,
    #                 1,
    #                 inputViewAngleImg,
    #                 0,
    #             )

    #             fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_LANDSAT8)
    #             fmaskConfig.setTOARefScaling(float(scaleFactor))
    #             fmaskConfig.setThermalInfo(thermalInfo)
    #             fmaskConfig.setAnglesInfo(anglesInfo)
    #             fmaskConfig.setKeepIntermediates(False)
    #             fmaskConfig.setVerbose(True)
    #             fmaskConfig.setTempDir(tmpBaseDIR)
    #             fmaskConfig.setMinCloudSize(minCloudSize)
    #             fmaskConfig.setEqn17CloudProbThresh(
    #                 fmask.config.FmaskConfig.Eqn17CloudProbThresh
    #             )
    #             fmaskConfig.setEqn20NirSnowThresh(
    #                 fmask.config.FmaskConfig.Eqn20NirSnowThresh
    #             )
    #             fmaskConfig.setEqn20GreenSnowThresh(
    #                 fmask.config.FmaskConfig.Eqn20GreenSnowThresh
    #             )

    #             # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
    #             toaImgInfo = rios.fileinfo.ImageInfo(inputReflImage)
    #             fmaskConfig.setCloudBufferSize(
    #                 int(cloudBufferDistance / toaImgInfo.xRes)
    #             )
    #             fmaskConfig.setShadowBufferSize(
    #                 int(shadowBufferDistance / toaImgInfo.xRes)
    #             )

    #             fmask.fmask.doFmask(fmaskFilenames, fmaskConfig)

    #             rsgislib.imagecalc.image_math(
    #                 tmpFMaskOut,
    #                 outputImage,
    #                 "(b1==2)?1:(b1==3)?2:0",
    #                 outFormat,
    #                 rsgislib.TYPE_8UINT,
    #             )

    #         elif cloud_msk_methods == "LSMSK":
    #             if (self.bandQAFile == "") or (not os.path.exists(self.bandQAFile)):
    #                 raise ARCSIException(
    #                     "The QA band is not present - cannot use this for cloud masking."
    #                 )

    #             bqa_img_file = self.bandQAFile
    #             if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
    #                 bqa_img_file, inputReflImage
    #             ):
    #                 bqa_img_file = os.path.join(tmpBaseDIR, tmpBaseName + "_BQA.kea")
    #                 rsgislib.imageutils.resample_img_to_match(
    #                     inputReflImage,
    #                     self.bandQAFile,
    #                     bqa_img_file,
    #                     "KEA",
    #                     rsgislib.INTERP_NEAREST_NEIGHBOUR,
    #                     rsgislib.TYPE_16UINT,
    #                     no_data_val=0,
    #                     multicore=False,
    #                 )

    #             if self.collection_num == 1:
    #                 exp = (
    #                     "(b1==2800)||(b1==2804)||(b1==2808)||(b1==2812)||(b1==6896)||(b1==6900)||(b1==6904)||(b1==6908)?1:"
    #                     "(b1==2976)||(b1==2980)||(b1==2984)||(b1==2988)||(b1==3008)||(b1==3012)||(b1==3016)||(b1==3020)||"
    #                     "(b1==7072)||(b1==7076)||(b1==7080)||(b1==7084)||(b1==7104)||(b1==7108)||(b1==7112)||(b1==7116)?2:0"
    #                 )
    #                 rsgislib.imagecalc.image_math(
    #                     bqa_img_file, outputImage, exp, outFormat, rsgislib.TYPE_8UINT
    #                 )
    #             elif self.collection_num == 2:
    #                 import rsgislib.imagecalibration.sensorlvl2data

    #                 c2_bqa_ind_img_file = os.path.join(
    #                     tmpBaseDIR, tmpBaseName + "c2_qa_ind_bands.kea"
    #                 )
    #                 rsgislib.imagecalibration.sensorlvl2data.parse_landsat_c2_qa_pixel_img(
    #                     bqa_img_file, c2_bqa_ind_img_file, gdalformat="KEA"
    #                 )
    #                 band_defns = list()
    #                 band_defns.append(
    #                     rsgislib.imagecalc.BandDefn(
    #                         "DilatedCloud", c2_bqa_ind_img_file, 2
    #                     )
    #                 )
    #                 band_defns.append(
    #                     rsgislib.imagecalc.BandDefn("Cloud", c2_bqa_ind_img_file, 4)
    #                 )
    #                 band_defns.append(
    #                     rsgislib.imagecalc.BandDefn(
    #                         "CloudShadow", c2_bqa_ind_img_file, 5
    #                     )
    #                 )
    #                 rsgislib.imagecalc.band_math(
    #                     outputImage,
    #                     "(DilatedCloud == 1)||(Cloud == 1)?1:(CloudShadow == 1)?2:0",
    #                     "KEA",
    #                     rsgislib.TYPE_8UINT,
    #                     band_defns,
    #                 )
    #             else:
    #                 raise ARCSIException(
    #                     "Can only read Collection 1 and 2 cloud masks."
    #                 )

    #         else:
    #             raise ARCSIException(
    #                 "Landsat only has FMASK and LSMSK cloud masking options; option provided is unknown."
    #             )

    #         if outFormat == "KEA":
    #             rsgislib.rastergis.pop_rat_img_stats(outputImage, True, True)
    #             ratDataset = gdal.Open(outputImage, gdal.GA_Update)
    #             red = rat.readColumn(ratDataset, "Red")
    #             green = rat.readColumn(ratDataset, "Green")
    #             blue = rat.readColumn(ratDataset, "Blue")
    #             ClassName = numpy.empty_like(red, dtype=numpy.dtype("a255"))

    #             red[0] = 0
    #             green[0] = 0
    #             blue[0] = 0

    #             if (red.shape[0] == 2) or (red.shape[0] == 3):
    #                 red[1] = 0
    #                 green[1] = 0
    #                 blue[1] = 255
    #                 ClassName[1] = "Clouds"

    #                 if red.shape[0] == 3:
    #                     red[2] = 0
    #                     green[2] = 255
    #                     blue[2] = 255
    #                     ClassName[2] = "Shadows"

    #             rat.writeColumn(ratDataset, "Red", red)
    #             rat.writeColumn(ratDataset, "Green", green)
    #             rat.writeColumn(ratDataset, "Blue", blue)
    #             rat.writeColumn(ratDataset, "ClassName", ClassName)
    #             ratDataset = None
    #         rsgislib.imageutils.copy_proj_from_img(outputImage, inputReflImage)

    #         if not self.debugMode:
    #             if not tmpDIRExisted:
    #                 shutil.rmtree(tmpBaseDIR, ignore_errors=True)

    #         return outputImage, None
    #     except Exception as e:
    #         raise e

    # def createCloudMaskDataArray(self, inImgDataArr):
    #     return inImgDataArr

    # def defineDarkShadowImageBand(self):
    #     return 5

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 7

    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        sixsCoeffs = numpy.zeros((3, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith #TODO
        s.geometry.solar_a = self.solarAzimuth #TODO
        s.geometry.view_z = self.sensorZenith #TODO
        s.geometry.view_a = self.sensorAzimuth #TODO
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = (
            float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
        )
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Blue
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["BAND5"].wvLenMin6S,
            self.specBandInfo["BAND5"].wvLenMax6S,
            self.specBandInfo["BAND5"].respFunc6S,
        )
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[0, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[0, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[0, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[0, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[0, 5] = float(s.outputs.values["environmental_irradiance"])

        # Green
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["BAND6"].wvLenMin6S,
            self.specBandInfo["BAND6"].wvLenMax6S,
            self.specBandInfo["BAND6"].respFunc6S,
        )
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[1, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[1, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[1, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[1, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[1, 5] = float(s.outputs.values["environmental_irradiance"])

        # Red
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["BAND7"].wvLenMin6S,
            self.specBandInfo["BAND7"].wvLenMax6S,
            self.specBandInfo["BAND7"].respFunc6S,
        )
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[2, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[2, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[2, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[2, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[2, 5] = float(s.outputs.values["environmental_irradiance"])

        # NIR BAND8
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["BAND8"].wvLenMin6S,
            self.specBandInfo["BAND8"].wvLenMax6S,
            self.specBandInfo["BAND8"].respFunc6S,
        )
        s.run()
        sixsCoeffs[3, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[3, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[3, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[3, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[3, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[3, 5] = float(s.outputs.values["environmental_irradiance"])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(
        self,
        inputRadImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        surfaceAltitude,
        aotVal,
        useBRDF,
        scaleFactor,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(
            aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
        )

        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=1,
                aX=float(sixsCoeffs[0, 0]),
                bX=float(sixsCoeffs[0, 1]),
                cX=float(sixsCoeffs[0, 2]),
                DirIrr=float(sixsCoeffs[0, 3]),
                DifIrr=float(sixsCoeffs[0, 4]),
                EnvIrr=float(sixsCoeffs[0, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=2,
                aX=float(sixsCoeffs[1, 0]),
                bX=float(sixsCoeffs[1, 1]),
                cX=float(sixsCoeffs[1, 2]),
                DirIrr=float(sixsCoeffs[1, 3]),
                DifIrr=float(sixsCoeffs[1, 4]),
                EnvIrr=float(sixsCoeffs[1, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=3,
                aX=float(sixsCoeffs[2, 0]),
                bX=float(sixsCoeffs[2, 1]),
                cX=float(sixsCoeffs[2, 2]),
                DirIrr=float(sixsCoeffs[2, 3]),
                DifIrr=float(sixsCoeffs[2, 4]),
                EnvIrr=float(sixsCoeffs[2, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=4,
                aX=float(sixsCoeffs[3, 0]),
                bX=float(sixsCoeffs[3, 1]),
                cX=float(sixsCoeffs[3, 2]),
                DirIrr=float(sixsCoeffs[3, 3]),
                DifIrr=float(sixsCoeffs[3, 4]),
                EnvIrr=float(sixsCoeffs[3, 5]),
            )
        )

        rsgislib.imagecalibration.apply_6s_coeff_single_param(
            inputRadImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            imgBandCoeffs,
        )
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotVal,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        scaleFactor,
        elevCoeffs=None,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                aotVal,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
            )
            print("LUT has been built.")

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=1,
                        aX=float(sixsCoeffs[0, 0]),
                        bX=float(sixsCoeffs[0, 1]),
                        cX=float(sixsCoeffs[0, 2]),
                        DirIrr=float(sixsCoeffs[0, 3]),
                        DifIrr=float(sixsCoeffs[0, 4]),
                        EnvIrr=float(sixsCoeffs[0, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=2,
                        aX=float(sixsCoeffs[1, 0]),
                        bX=float(sixsCoeffs[1, 1]),
                        cX=float(sixsCoeffs[1, 2]),
                        DirIrr=float(sixsCoeffs[1, 3]),
                        DifIrr=float(sixsCoeffs[1, 4]),
                        EnvIrr=float(sixsCoeffs[1, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=3,
                        aX=float(sixsCoeffs[2, 0]),
                        bX=float(sixsCoeffs[2, 1]),
                        cX=float(sixsCoeffs[2, 2]),
                        DirIrr=float(sixsCoeffs[2, 3]),
                        DifIrr=float(sixsCoeffs[2, 4]),
                        EnvIrr=float(sixsCoeffs[2, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=4,
                        aX=float(sixsCoeffs[3, 0]),
                        bX=float(sixsCoeffs[3, 1]),
                        cX=float(sixsCoeffs[3, 2]),
                        DirIrr=float(sixsCoeffs[3, 3]),
                        DifIrr=float(sixsCoeffs[3, 4]),
                        EnvIrr=float(sixsCoeffs[3, 5]),
                    )
                )

                elevCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=imgBandCoeffs
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_lut_param(
            inputRadImage,
            inputDEMFile,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevCoeffs,
        )
        return outputImage, elevCoeffs

    def convertImageToSurfaceReflAOTDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        inputAOTImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        aotMin,
        aotMax,
        scaleFactor,
        elevAOTCoeffs=None,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
                aotMin,
                aotMax,
            )

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=1,
                            aX=float(sixsCoeffs[0, 0]),
                            bX=float(sixsCoeffs[0, 1]),
                            cX=float(sixsCoeffs[0, 2]),
                            DirIrr=float(sixsCoeffs[0, 3]),
                            DifIrr=float(sixsCoeffs[0, 4]),
                            EnvIrr=float(sixsCoeffs[0, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=2,
                            aX=float(sixsCoeffs[1, 0]),
                            bX=float(sixsCoeffs[1, 1]),
                            cX=float(sixsCoeffs[1, 2]),
                            DirIrr=float(sixsCoeffs[1, 3]),
                            DifIrr=float(sixsCoeffs[1, 4]),
                            EnvIrr=float(sixsCoeffs[1, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=3,
                            aX=float(sixsCoeffs[2, 0]),
                            bX=float(sixsCoeffs[2, 1]),
                            cX=float(sixsCoeffs[2, 2]),
                            DirIrr=float(sixsCoeffs[2, 3]),
                            DifIrr=float(sixsCoeffs[2, 4]),
                            EnvIrr=float(sixsCoeffs[2, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=4,
                            aX=float(sixsCoeffs[3, 0]),
                            bX=float(sixsCoeffs[3, 1]),
                            cX=float(sixsCoeffs[3, 2]),
                            DirIrr=float(sixsCoeffs[3, 3]),
                            DifIrr=float(sixsCoeffs[3, 4]),
                            EnvIrr=float(sixsCoeffs[3, 5]),
                        )
                    )
                    aot6SCoeffsOut.append(
                        rsgislib.imagecalibration.AOTLUTFeat(
                            AOT=float(aotVal), Coeffs=imgBandCoeffs
                        )
                    )
                elevAOTCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=aot6SCoeffsOut
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_aot_lut_param(
            inputRadImage,
            inputDEMFile,
            inputAOTImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevAOTCoeffs,
        )

        return outputImage, elevAOTCoeffs

    # def run6SToOptimiseAODValue(
    #     self,
    #     aotVal,
    #     radBlueVal,
    #     predBlueVal,
    #     aeroProfile,
    #     atmosProfile,
    #     grdRefl,
    #     surfaceAltitude,
    # ):

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        raise ARCSIException("Not Implemented")

    def estimateImageToAODUsingDDV(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
    ):
        raise ARCSIException("Not Implemented")

    def estimateImageToAODUsingDOS(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
        globalDOS,
        simpleDOS,
        dosOutRefl,
    ):
        raise ARCSIException("Not Implemented")

    def estimateSingleAOTFromDOS(
        self,
        radianceImage,
        toaImage,
        inputDEMFile,
        tmpPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        minAOT,
        maxAOT,
        dosOutRefl,
    ):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(
                radianceImage,
                toaImage,
                inputDEMFile,
                tmpPath,
                outputName,
                outFormat,
                aeroProfile,
                atmosProfile,
                grdRefl,
                minAOT,
                maxAOT,
                dosOutRefl,
                2,
            )
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("Blue")
        dataset.GetRasterBand(2).SetDescription("Green")
        dataset.GetRasterBand(3).SetDescription("Red")
        dataset.GetRasterBand(4).SetDescription("NIR")
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")
