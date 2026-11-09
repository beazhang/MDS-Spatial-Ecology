"""
Segmentation-mask helpers.

Library form of the helper used by the Baysor / 10x multimodal megakaryocyte
resegmentation. That resegmentation runs upstream of the figures, so no figure notebook
in this repository imports this module.
"""
import geopandas as gpd
import numpy as np
import shapely


def contours_to_polygons(
    x,
    y,
    ids,
) -> gpd.GeoDataFrame:

    x = np.asarray(x)
    y = np.asarray(y)
    ids = np.asarray(ids)
    
    splits = np.where(ids[:-1] != ids[1:])[0] + 1
    geometry_offset = np.hstack([0, splits, len(ids)])
    part_offset = np.arange(len(np.unique(ids)) + 1)
    polygons = shapely.from_ragged_array(
        shapely.GeometryType.POLYGON,
        coords=np.stack([x, y]).T.copy(order='C'),
        offsets=(geometry_offset, part_offset),
    )

    indices = np.sort(np.unique(ids, return_index=True)[1])
    return gpd.GeoDataFrame(geometry=polygons, index=ids[indices])
