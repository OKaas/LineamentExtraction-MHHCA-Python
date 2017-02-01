# lineament-extraction-python
A new algorithm for spatial clustering of the line segments as a tool for lineament extraction.

The lineament is a linear feature describing discontinuity in a landscape. The lineament extraction is not an easy problem. Recently, an automatic approach based on multi-hillshade hierarchic clustering has been developed. An important part of this approach is the spatial line segment clustering. This paper presents a new algorithm for this clustering, based on a facility location algorithm. The results of implementation show acceleration in comparison with its predecessor.

### Requirements: 
* python3 - check [Python](https://www.python.org/downloads/)

### Usage:
* ``$python3 __main__.py [-h] <input file path> <border X> <border Y> <output file path>``
    
###### Where:
* ``<input file path>``  - path to input file
*    ``<border X>``         - width of buffer zone
*    ``<border Y>``         - height of buffer zone
*    ``<border azimuth>``   - azimuth threshold
*    ``<output file path>`` - path to output file
    
###### Example:
*    ``$python3 __main__.py C:\Data\input.txt 150 200 20 C:\Data\result.txt``

### Input data:
* [MHHCA-Python-Input](https://github.com/OKaas/LineamentExtraction-MHHCA-Input)