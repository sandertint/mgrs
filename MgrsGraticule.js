/**
 * @module ol/layer/Graticule
 */
import Feature from 'ol/Feature.js';
import GeometryLayout from 'ol/geom/GeometryLayout.js';
import MultiLineString from 'ol/geom/MultiLineString.js';
import Point from 'ol/geom/Point.js';
import {
  applyTransform,
  containsCoordinate,
  containsExtent,
  containsXY,
  equals,
  getCenter,
  getHeight,
  getIntersection,
  getWidth,
  isEmpty,
} from 'ol/extent.js';
import {clamp} from 'ol/math';
import {
  equivalent as equivalentProjection,
} from 'ol/proj';
import {getVectorContext} from 'ol/render';

import proj4 from "proj4";
import { Graticule } from 'ol';
import Style from 'ol/style/Style';
import Stroke from 'ol/style/Stroke';
import Fill from 'ol/style/Fill';
import Text from 'ol/style/Text';

const DEFAULT_STROKE_STYLE_BOLD = new Stroke({
  color: 'rgba(0,0,0,0.2)',
  width: 3
});

class MgrsGraticule extends Graticule {
  constructor(opt_options) {
    const options = opt_options ? opt_options : {};

    const zones = "zones" in options ? options.zones : ["34V", "35V"];
    proj4.defs(
      zones.map((zone) => {
        const zoneNr = zone.replace(/\D+/g, "");
        return [
        `MGRS:${zoneNr}`, `+proj=utm +zone=${zoneNr} +datum=WGS84 +units=m +no_defs`
        ];
      })
    );

    options.showLabels = true;

    super(options);

    this.zones = zones;

    this.gridZoneDesignatorLabels_ = [];

    this.lineStyleBold_ = new Style({
      stroke: DEFAULT_STROKE_STYLE_BOLD,
    });

    this.zoneLabelStyleBase_ = new Style({
      text: new Text({
              font: '40px Calibri,sans-serif',
              textAlign: 'center',
              fill: new Fill({
                color: 'rgba(0,0,0,0.2)',
              }),
              stroke: new Stroke({
                color: 'rgba(255,255,255,0)',
                width: 0,
              }),
            }),
    });

    this.gridZoneLabelStyle_ = function (feature) {
      const label = feature.get('zone_label');
      this.zoneLabelStyleBase_.getText().setText(label);
      return this.zoneLabelStyleBase_;
    }.bind(this);
  }

  loaderFunction(extent, resolution, projection) {
    this.loadedExtent_ = extent;
    const source = this.getSource();

    // only consider the intersection between our own extent & the requested one
    const layerExtent = this.getExtent() || [
      -Infinity,
      -Infinity,
      Infinity,
      Infinity,
    ];
    const renderExtent = getIntersection(layerExtent, extent);

    if (this.renderedExtent_ && equals(this.renderedExtent_, renderExtent)) {
      return;
    }
    this.renderedExtent_ = renderExtent;

    // bail out if nothing to render
    if (isEmpty(renderExtent)) {
      return;
    }

    // update projection info
    const center = getCenter(renderExtent);
    const squaredTolerance = (resolution * resolution) / 4;

    const updateProjectionInfo =
      !this.projection_ || !equivalentProjection(this.projection_, projection);

    if (updateProjectionInfo) {
      this.updateProjectionInfo_(projection);
    }

    this.createGraticule_(renderExtent, center, resolution, squaredTolerance);

    // first make sure we have enough features in the pool
    let featureCount = this.meridians_.length + this.parallels_.length;
    if (this.meridiansLabels_) {
      featureCount += this.meridians_.length;
    }
    if (this.parallelsLabels_) {
      featureCount += this.parallels_.length;
    }

    if (this.gridZoneDesignatorLabels_) {
      featureCount += this.gridZoneDesignatorLabels_.length;
    }

    let feature;
    while (featureCount > this.featurePool_.length) {
      feature = new Feature();
      this.featurePool_.push(feature);
    }

    const featuresColl = source.getFeaturesCollection();
    featuresColl.clear();
    let poolIndex = 0;

    // add features for the lines & labels
    let i, l;
    for (i = 0, l = this.meridians_.length; i < l; ++i) {
      feature = this.featurePool_[poolIndex++];
      feature.setGeometry(this.meridians_[i]);
      feature.setStyle(!this.meridians_[i].get("isGridZone") ? this.lineStyle_ : this.lineStyleBold_);
      featuresColl.push(feature);
    }
    for (i = 0, l = this.parallels_.length; i < l; ++i) {
      feature = this.featurePool_[poolIndex++];
      feature.setGeometry(this.parallels_[i]);
      feature.setStyle(!this.parallels_[i].get("isGridZone") ? this.lineStyle_ : this.lineStyleBold_);
      featuresColl.push(feature);
    }
  }

  createGraticule_(extent, center, resolution) {
    const interval = this.getInterval_(resolution);

    if (interval == -1) {
      this.meridians_.length = 0;
      this.parallels_.length = 0;console.log(feature);
      if (this.meridiansLabels_) {
        this.meridiansLabels_.length = 0;
      }
      if (this.parallelsLabels_) {
        this.parallelsLabels_.length = 0;
      }
      return;
    }

    let wrapX = false;
    const projectionExtent = this.projection_.getExtent();
    const worldWidth = getWidth(projectionExtent);
    if (
      this.getSource().getWrapX() &&
      this.projection_.canWrapX() &&
      !containsExtent(projectionExtent, extent)
    ) {
      if (getWidth(extent) >= worldWidth) {
        extent[0] = projectionExtent[0];
        extent[2] = projectionExtent[2];
      } else {
        wrapX = true;
      }
    }

    // Constrain the center to fit into the extent available to the graticule

    const validCenterP = [
      clamp(center[0], this.minX_, this.maxX_),
      clamp(center[1], this.minY_, this.maxY_),
    ];

    // Transform the center to lon lat
    // Some projections may have a void area at the poles
    // so replace any NaN latitudes with the min or max value closest to a pole

    const centerLonLat = this.toLonLatTransform_(validCenterP);
    if (isNaN(centerLonLat[1])) {
      centerLonLat[1] =
        Math.abs(this.maxLat_) >= Math.abs(this.minLat_)
          ? this.maxLat_
          : this.minLat_;
    }
    let centerLon = clamp(centerLonLat[0], this.minLon_, this.maxLon_);
    let centerLat = clamp(centerLonLat[1], this.minLat_, this.maxLat_);

    // Limit the extent to fit into the extent available to the graticule

    let validExtentP = extent;
    if (!wrapX) {
      validExtentP = [
        clamp(extent[0], this.minX_, this.maxX_),
        clamp(extent[1], this.minY_, this.maxY_),
        clamp(extent[2], this.minX_, this.maxX_),
        clamp(extent[3], this.minY_, this.maxY_),
      ];
    }

    // Transform the extent to get the lon lat ranges for the edges of the extent

    const validExtent = applyTransform(
      validExtentP,
      this.toLonLatTransform_,
      undefined,
      8
    );

    let maxLat = validExtent[3];
    let maxLon = validExtent[2];
    let minLat = validExtent[1];
    let minLon = validExtent[0];

    if (!wrapX) {
      // Check if extremities of the world extent lie inside the extent
      // (for example the pole in a polar projection)
      // and extend the extent as appropriate

      if (containsCoordinate(validExtentP, this.bottomLeft_)) {
        minLon = this.minLon_;
        minLat = this.minLat_;
      }
      if (containsCoordinate(validExtentP, this.bottomRight_)) {
        maxLon = this.maxLon_;
        minLat = this.minLat_;
      }
      if (containsCoordinate(validExtentP, this.topLeft_)) {
        minLon = this.minLon_;
        maxLat = this.maxLat_;
      }
      if (containsCoordinate(validExtentP, this.topRight_)) {
        maxLon = this.maxLon_;
        maxLat = this.maxLat_;
      }

      // The transformed center may also extend the lon lat ranges used for rendering

      maxLat = clamp(maxLat, centerLat, this.maxLat_);
      maxLon = clamp(maxLon, centerLon, this.maxLon_);
      minLat = clamp(minLat, this.minLat_, centerLat);
      minLon = clamp(minLon, this.minLon_, centerLon);
    }

    let gridSize = 100000;

    if (interval === 0.5) {
      gridSize = 50000;
    } else if (interval <= 0.2 && interval > 0.05) {
      gridSize = 10000;
    } else if (interval <= 0.05 && interval > 0.01) {
      gridSize = 5000;
    } else if (interval <= 0.01) {
      gridSize = 1000;
    }

    let meridianCount = 0;
    let parallelCount = 0;
    let gridZoneLabelsCount = 0;
    this.zones.forEach((zone) => {
      const { numAddedMeridians, numAddedParallels, numAddedGridZoneLabels } =
        this.getZoneGraticule_(
          zone,
          [minLon, minLat, maxLon, maxLat],
          gridSize,
          meridianCount,
          parallelCount,
          gridZoneLabelsCount
        );
      meridianCount += numAddedMeridians;
      parallelCount += numAddedParallels;
      gridZoneLabelsCount += numAddedGridZoneLabels;
    });

    this.meridians_.length = meridianCount;
    this.meridiansLabels_.length = meridianCount;
    this.parallels_.length = parallelCount;
    this.parallelsLabels_.length = parallelCount;
    this.gridZoneDesignatorLabels_.length = gridZoneLabelsCount;
  }

  getZoneGraticule_(zoneCode, extent, gridSize, meridianCount, parallelCount, gridZoneLabelsCount) {

    const zoneLatCode = zoneCode.slice(-1);
    const zoneLonCode = zoneCode.slice(0, 2);
    const maxLonDeg = parseInt(zoneLonCode) * 6 - 180;
    
    const minLatDeg = ["N", "P", "Q", "R", "S", "T", "U", "V", "W"].indexOf(zoneLatCode) * 8;
    const [minLonDeg, maxLatDeg] = [maxLonDeg - 6, minLatDeg + 8];

    const zoneSrs = `MGRS:${zoneLonCode}`;

    const view_extent_m = [
      ...proj4("EPSG:4326", zoneSrs, [extent[0], extent[1]]),
      ...proj4("EPSG:4326", zoneSrs, [extent[2], extent[3]])
    ];

    const x_min = (view_extent_m[0] - (view_extent_m[0] % gridSize)) - gridSize;
    const x_max = (view_extent_m[2] - (view_extent_m[2] % gridSize)) + gridSize;
 
    const y_min = (view_extent_m[1] - (view_extent_m[1] % gridSize)) - gridSize;
    const y_max = (view_extent_m[3] - (view_extent_m[3] % gridSize)) + gridSize;
    
    const num_x = (x_max - x_min) / gridSize;
    const num_y = (y_max - y_min) / gridSize;
 
    const minEpsg3857 = proj4("EPSG:4326", "EPSG:3857", [minLonDeg, minLatDeg]);
    const maxEpsg3857 = proj4("EPSG:4326", "EPSG:3857", [maxLonDeg, maxLatDeg]);

    const extentEpsg3857 = [...minEpsg3857, ...maxEpsg3857];
    
    let i = 0;
    for (; i <= num_x; i++) {
      const flatCoordinates = [];
      for (let j=0; j <= num_y; j++) {
        const coords = proj4(zoneSrs, "EPSG:3857", [
          x_min + (gridSize * i),
          y_min + (gridSize * j)
        ]);

        flatCoordinates.push([
          clamp(coords[0], minEpsg3857[0], maxEpsg3857[0]),
          clamp(coords[1], minEpsg3857[1], maxEpsg3857[1])
        ]);
      }

      let lineString = this.meridians_[i + meridianCount];
      if (!lineString) {
        lineString = new MultiLineString([flatCoordinates], GeometryLayout.XY);
        this.meridians_.push(lineString);
      } else {
        lineString.setCoordinates([flatCoordinates]);
        lineString.changed();
      }

      lineString.set("isGridZone", (x_min + (gridSize * i)) % 100000 === 0);

      if (this.meridiansLabels_) {
        if (i + meridianCount in this.meridiansLabels_) {
          this.meridiansLabels_[i + meridianCount].text = x_min + (gridSize * i);
        } else {
          this.meridiansLabels_[i + meridianCount] = {
            geom: new Point([]),
            text: x_min + (gridSize * i),
          };
        }
      }
    }

    // Create parallels
    let ii=0;
    for (; ii <= num_y; ii++) {
      const flatCoordinates = [];
      for (let jj=0; jj <= num_x; jj++) {
        const [x, y] = [x_min + (gridSize * jj), y_min + (gridSize * ii)];
        const coords1 = proj4(zoneSrs, "EPSG:3857", [x, y]);

        const x3857 = clamp(coords1[0], minEpsg3857[0], maxEpsg3857[0]);
        const y3857 = clamp(coords1[1], minEpsg3857[1], maxEpsg3857[1]);
        
        const prevCoord = flatCoordinates.length > 0 && flatCoordinates[flatCoordinates.length - 1];
        const isSameX = prevCoord && x3857.toPrecision(8) === prevCoord[0].toPrecision(8);

        if (prevCoord && isSameX) {
          if (flatCoordinates.length === 1) {
            flatCoordinates[0] = [x3857, y3857];
          }
        } else {
          flatCoordinates.push([x3857, y3857]);
        }
      }

      let lineString = this.parallels_[ii + parallelCount];
      if (!lineString) {
        lineString = new MultiLineString([flatCoordinates], GeometryLayout.XY);
        this.parallels_.push(lineString);
      } else {
        lineString.setCoordinates([flatCoordinates]);
        lineString.changed();
      }

      lineString.set(
        "isGridZone",
        (y_min + (gridSize * ii)) % 100000 === 0
      );

        if (this.parallelsLabels_) {
          let substrSize = 2;
          if (gridSize === 10000) {
            substrSize = 2;
          } else if (gridSize === 5000) {
            substrSize = 3;
          } else if (gridSize === 1000) {
            substrSize = 4;
          }
          const text = `${(y_min + (gridSize * ii))}`.substring(-substrSize);
          
          if (ii + parallelCount in this.parallelsLabels_) {
            this.parallelsLabels_[ii + parallelCount].text = text;
          } else {
            this.parallelsLabels_[ii + parallelCount] = {
              geom: new Point([]),
              text,
            };
          }
        }
    }

    const numAddedGridZoneLabels = this.getZoneCodeLabels_(
      view_extent_m, zoneLonCode, zoneSrs, extentEpsg3857, gridZoneLabelsCount);

    return { numAddedMeridians: i, numAddedParallels: ii, numAddedGridZoneLabels };
  }

  getZoneCodeLabels_(view_extent_m, zoneLonCode, zoneSrs, extentEpsg3857, gridZoneLabelsCount) {
    const gridSize = 100000;
    const x_min = (view_extent_m[0] - (view_extent_m[0] % gridSize));
    const x_max = (view_extent_m[2] - (view_extent_m[2] % gridSize));

    const y_min = (view_extent_m[1] - (view_extent_m[1] % gridSize));
    const y_max = (view_extent_m[3] - (view_extent_m[3] % gridSize));

    const num_x = (x_max - x_min) / gridSize;
    const num_y = (y_max - y_min) / gridSize;

    let idx = 0;
    for (let i = 0; i <= num_x; i++) {
      const curr_delta_x = x_min + i * gridSize;
      let x = curr_delta_x + (gridSize / 2);
      const startsOutsideExtentLeft = curr_delta_x < view_extent_m[0];
      const endsOutsideExtentRight = (curr_delta_x + gridSize) > view_extent_m[2];
      if (startsOutsideExtentLeft && endsOutsideExtentRight) {
        x = (view_extent_m[0] + view_extent_m[2]) / 2;
      } else if (startsOutsideExtentLeft) {
        x = ((x_min + gridSize) + view_extent_m[0]) / 2;
      } else if (endsOutsideExtentRight) {
        x = (curr_delta_x + view_extent_m[2]) / 2;
      }

      for (let j = 0; j <= num_y; j++) {
        const curr_delta_y = y_min + j * gridSize;
        let y = curr_delta_y + (gridSize / 2);
        const startsOutsideExtentBottom = curr_delta_y < view_extent_m[1];
        const endsOutsideExtentTop = curr_delta_y + gridSize > view_extent_m[3];
        if (startsOutsideExtentBottom && endsOutsideExtentTop) {
          y = (view_extent_m[1] + view_extent_m[3]) / 2; 
        } else if (startsOutsideExtentBottom) {
          y = ((y_min + gridSize) + view_extent_m[1]) / 2;
        } else if (endsOutsideExtentTop) {
          y = (curr_delta_y + view_extent_m[3]) / 2;
        }

        const coords2 = proj4(zoneSrs, "EPSG:3857", [x, y]);

        if (!containsXY(extentEpsg3857, ...coords2)) {
          continue;
        }

        const x3857 = clamp(coords2[0], extentEpsg3857[0], extentEpsg3857[2]);
        const y3857 = clamp(coords2[1], extentEpsg3857[1], extentEpsg3857[3]);

        const zoneCode = this.getZoneCode_(x, y, +zoneLonCode);

        const offsetIdx = idx + gridZoneLabelsCount;
        const gridZoneLabel = this.gridZoneDesignatorLabels_[offsetIdx];

        if (typeof gridZoneLabel !== "undefined") {
          this.gridZoneDesignatorLabels_[offsetIdx].text = zoneCode;
          this.gridZoneDesignatorLabels_[offsetIdx].geom.setCoordinates([x3857, y3857]);
        } else {
          this.gridZoneDesignatorLabels_[offsetIdx] = {
            geom: new Point([x3857, y3857]),
            text: zoneCode
          }
        }
        idx++;
      }
    }
    return idx;
  }

  getZoneCode_ (x, y, zoneNr) {
    const chars = "ABCDEFGHJKLMNPQRSTUVWXYZ";
    const charsSectored = ["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"];
    const latStartsWithF = zoneNr % 2 === 0;
    const lonType = zoneNr % 3;
    const lonTypeSector = charsSectored[lonType - 1];
    return lonTypeSector[+(`${x}`[0]) - 1] + chars[+(`${y}`[1]) + (latStartsWithF ? 5 : 0)];
  }

  drawLabels_(event) {
    const rotation = event.frameState.viewState.rotation;
    const extent = event.frameState.extent;
    const rotationCenter = getCenter(extent);
    let rotationExtent = extent;
    if (rotation) {
      const width = getWidth(extent);
      const height = getHeight(extent);
      const cr = Math.abs(Math.cos(rotation));
      const sr = Math.abs(Math.sin(rotation));
      const unrotatedWidth = (sr * height - cr * width) / (sr * sr - cr * cr);
      const unrotatedHeight = (sr * width - cr * height) / (sr * sr - cr * cr);
      rotationExtent = [
        rotationCenter[0] - unrotatedWidth / 2,
        rotationCenter[1] - unrotatedHeight / 2,
        rotationCenter[0] + unrotatedWidth / 2,
        rotationCenter[1] + unrotatedHeight / 2,
      ];
    }

    let startWorld = 0;
    let endWorld = 0;
    let labelsAtStart = this.latLabelPosition_ < 0.5;
    const projectionExtent = this.projection_.getExtent();
    const worldWidth = getWidth(projectionExtent);
    if (
      this.getSource().getWrapX() &&
      this.projection_.canWrapX() &&
      !containsExtent(projectionExtent, extent)
    ) {
      startWorld = Math.floor((extent[0] - projectionExtent[0]) / worldWidth);
      endWorld = Math.ceil((extent[2] - projectionExtent[2]) / worldWidth);
      const inverted = Math.abs(rotation) > Math.PI / 2;
      labelsAtStart = labelsAtStart !== inverted;
    }
    const vectorContext = getVectorContext(event);

    for (let world = startWorld; world <= endWorld; ++world) {
      let poolIndex = this.meridians_.length + this.parallels_.length;
      let feature, index, l, textPoint;

      if (this.meridiansLabels_) {
        for (index = 0, l = this.meridiansLabels_.length; index < l; ++index) {
          const lineString = this.meridians_[index];
          if (!rotation && world === 0) {
            textPoint = this.getMeridianPoint_(lineString, extent, index);
          } else {
            const clone = lineString.clone();
            clone.translate(world * worldWidth, 0);
            clone.rotate(-rotation, rotationCenter);
            textPoint = this.getMeridianPoint_(clone, rotationExtent, index);
            textPoint.rotate(rotation, rotationCenter);
          }
          feature = this.featurePool_[poolIndex++];
          feature.setGeometry(textPoint);
          feature.set('graticule_label', this.meridiansLabels_[index].text);
          vectorContext.drawFeature(feature, this.lonLabelStyle_(feature));
        }
      }
      if (this.parallelsLabels_) {
        if (
          (world === startWorld && labelsAtStart) ||
          (world === endWorld && !labelsAtStart)
        ) {
          for (index = 0, l = this.parallels_.length; index < l; ++index) {
            const lineString = this.parallels_[index];
            if (!rotation && world === 0) {
              textPoint = this.getParallelPoint_(lineString, extent, index);
            } else {
              const clone = lineString.clone();
              clone.translate(world * worldWidth, 0);
              clone.rotate(-rotation, rotationCenter);
              textPoint = this.getParallelPoint_(clone, rotationExtent, index);
              textPoint.rotate(rotation, rotationCenter);
            }
            feature = this.featurePool_[poolIndex++];
            feature.setGeometry(textPoint);
            feature.set('graticule_label', this.parallelsLabels_[index].text);
            vectorContext.drawFeature(feature, this.latLabelStyle_(feature));
          }
        }
      }

      if (this.gridZoneDesignatorLabels_) {
        const gridZoneLabels = [...this.gridZoneDesignatorLabels_];
        for (index = 0, l = gridZoneLabels.length; index < l; ++index) {
          feature = this.featurePool_[poolIndex++];
          feature.setGeometry(gridZoneLabels[index].geom);
          feature.set('zone_label', gridZoneLabels[index].text);
          feature.setStyle(this.gridZoneLabelStyle_);
        
          vectorContext.drawFeature(feature, this.gridZoneLabelStyle_(feature));
        }
      }
    }
  }
}

export default MgrsGraticule;