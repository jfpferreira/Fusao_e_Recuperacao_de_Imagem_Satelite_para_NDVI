var L7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    Modis = ee.ImageCollection("MODIS/006/MOD09GQ"),
    osmwater = ee.FeatureCollection("users/jfpferreira/OSM_Water"),
    osmroads = ee.FeatureCollection("users/jfpferreira/OSM_Roads"),
    osmbuildings = ee.FeatureCollection("users/jfpferreira/OSM_Buildings");

var saocristovao = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '181621');
var santar = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '180911');
var benfica = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '140302');
var beringel = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '020503');

var ndviParams = {
    min: 0.20,
    max: 0.80,
    palette: ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#006837', '#004529']
}
var errorParams = {
    min: -0.2,
    max: 0.2,
    palette: ['#fc8d59', '#ffffbf', '#91cf60']
}

// Collection
var collection = ee.ImageCollection(ee.Image.constant(0));
// to Modis Imagery
var collectionModis = ee.ImageCollection(ee.Image.constant(0));
// Imagery to recover
var collectionR = ee.ImageCollection(ee.Image.constant(0));


var maskClouds = function(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = ee.Number(2).pow(3).int();
    var cloudsBitMask = ee.Number(2).pow(5).int();
    // Get the pixel QA band.
    var qa = image.select('pixel_qa');
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    return image.updateMask(mask)

};

// function to add NDVI and time
var addNDVI = function(img) {
    return img.addBands(img.normalizedDifference(['nir', 'red']).rename('NDVI'))
};

var addTime = function(img) {
    return img
        // System Time
        .addBands(img.metadata('system:time_start').rename('Time')).set('Time', img.get('system:time_start'))
};

var linearRegression = function(image, band, roi, windowBefore, windowAfter) {

    var time_img = image.get('system:time_start')

    // imagery start date
    var startWindowDate = ee.Date(ee.Number(time_img).subtract(windowBefore * (1000 * 60 * 60 * 24)));
    // imagery end date
    var endWindowDate = ee.Date(ee.Number(time_img).add(windowAfter * (1000 * 60 * 60 * 24)));

    // Only imgery that have less than 5% of cloud cover in that region
    var colLR = collection.filterDate(startWindowDate, endWindowDate)
        .filter(ee.Filter.notEquals('system:time_start', time_img))
        .filterMetadata('CLOUD_COVER_REGION', 'less_than', 5)

    var lr = colLR.select('Time', band)
        .reduce(ee.Reducer.linearFit());

    var a = lr.select('scale');
    var b = lr.select('offset');
    var predicted_lr = a.multiply(ee.Number(time_img)).add(b).rename(band)

    //Mean of the estimated pixels
    var mean_pixels = predicted_lr.select(band).reduceRegion(ee.Reducer.mean(), roi, 30).toImage([band]).rename('mean_lr')

    return image.select(band).addBands(image.metadata('system:time_start').rename('Time')).addBands(mean_pixels.rename('mean_lr'))
        .addBands(predicted_lr.select(band).rename('estimated')).unmask(predicted_lr)
        .set('LR_Size', colLR.size())

}

// Perform Global correction
var globalCorrection = function(image, roi, band, windowBefore, windowAfter) {

    // Original Image - T0
    // Cf pixels - outside the stripes
    var cf_t0 = image.reduceRegion(ee.Reducer.mean(), roi).toImage([band])

    //Closest image from original Image - T1
    var t1 = collection.map(function(img) {
            var time = ee.Number(image.get('system:time_start'));
            var dif = ee.Number(img.get('system:time_start')).subtract(time).abs()
            return img.set('Distance', dif)
        }).filter(ee.Filter.notEquals('system:time_start', ee.Number(image.get('system:time_start'))))
        .filterMetadata('CLOUD_COVER_REGION', 'equals', 0).sort('Distance', true).first()

    var shortTimeFrame = t1.select(band)

    //Closest Image - T1
    // Want Cg pixels - neighboor pixels inside the stripes (eq(-2))
    var shortTimeFrame_strips = t1.mask(image.select(band).unmask(-2).eq(-2))
    var cg_t1 = shortTimeFrame_strips.reduceRegion(ee.Reducer.mean(), roi).toImage([band])

    // Want Cf pixels - neighboor pixels outside the stripes (neq(-2))
    var shortTimeFrame_nonStrips = t1.mask(image.select(band).unmask(-2).neq(-2)).updateMask(shortTimeFrame)
    var cf_t1 = shortTimeFrame_nonStrips.reduceRegion(ee.Reducer.mean(), roi).toImage([band])

    // Delta - Difference between Cg from T1 and Cf
    var delta = cg_t1.subtract(cf_t1)

    // Linear Regression Image 
    var lr = linearRegression(image, band, roi, windowBefore, windowAfter)
    var cg_linear = lr.select('estimated').mask(image.select(band).unmask(-2).eq(-2)).reduceRegion(ee.Reducer.mean(), roi, 30).toImage(['estimated'])

    var pixels_toRecover_cg = image.select(band).unmask(lr.select(band).add(cf_t0).add(delta).subtract(cg_linear)).rename(band)

    return pixels_toRecover_cg.addBands(image.metadata('system:time_start').rename('Time'))
        .addBands(image.select(band).subtract(pixels_toRecover_cg.select(band)).rename('error_LST'))
        .copyProperties(lr, ['LR_Size'])

}

// To get mean for each cluster
var groupCluster = function(srcimg, dstimg, roi, band, name) {

    var group = ee.List(srcimg.select(band, 'cluster').reduceRegion({
        reducer: ee.Reducer.mean().group({
            groupField: 1,
            groupName: 'cluster',
        }),
        geometry: roi,
        scale: 30
    }).get('groups'))

    // Get the cluster keys and their corresponding values
    var keys = group.map(function(item) {
        return ee.Dictionary(item).get('cluster')
    })
    var vals = group.map(function(item) {
        return ee.Dictionary(item).get('mean')
    })

    return ee.Algorithms.If(ee.Algorithms.IsEqual(group, ee.List([])),
        dstimg.addBands(ee.Image(0).rename(name)),
        dstimg.select('cluster').remap(keys, vals).rename(name))

}

// Perform Adaptive correction
var adaptiveCorrection = function(image, roi, k, band, windowBefore, windowAfter) {

    // Make the training dataset.
    var training = image.select(band).sample(roi, 30);
    var clusterer = ee.Clusterer.wekaKMeans(k).train(training);

    //Strips Image Original - T0
    // Want Cf pixels - neighboor pixels outside the stripes
    // It is already masked
    var original = image.addBands(image.cluster(clusterer)) //add band 'cluster'

    //Closest image from original Image
    var t1 = collection.map(function(img) {
            var time = ee.Number(image.get('system:time_start'));
            var dif = ee.Number(img.get('system:time_start')).subtract(time).abs()
            return img.set('Distance', dif)
        }).filter(ee.Filter.notEquals('system:time_start', ee.Number(image.get('system:time_start'))))
        .filterMetadata('CLOUD_COVER_REGION', 'equals', 0).sort('Distance', true).first()

    var shortTimeFrame = t1.select(band).addBands(t1.cluster(clusterer))

    //Closest Image - T1
    // Want Cg pixels - neighboor pixels inside the stripes (eq(-2))
    var shortTimeFrame_strips = shortTimeFrame.mask(image.select(band).unmask(-2).eq(-2))
    var cg_t1 = shortTimeFrame_strips // with same training set

    // Want Cf pixels - neighboor pixels outside the stripes (neq(-2))
    var shortTimeFrame_nonStrips = shortTimeFrame.mask(image.select(band).unmask(-2).neq(-2))
    var cf_t1 = shortTimeFrame_nonStrips // with same training set

    // Linear Regression Image 
    var lr = linearRegression(image, band, roi, windowBefore, windowAfter)
        .addBands(linearRegression(image, band, roi, windowBefore, windowAfter).cluster(clusterer))

    var pixels_toRecover = lr.addBands(groupCluster(original, lr, roi, band, 'cf_t0'))
        .addBands(groupCluster(cg_t1, lr, roi, band, 'cg_t1'))
        .addBands(groupCluster(cf_t1, lr, roi, band, 'cf_t1'))
        .addBands(groupCluster(lr.mask(image.select(band).unmask(-2).eq(-2)), lr, roi, 'estimated', 'mean_cg_lr'))


    var delta = pixels_toRecover.addBands(pixels_toRecover.select('cg_t1').subtract(pixels_toRecover.select('cf_t1')).rename('delta'))
    var ndvi_cg = image.select(band).unmask(delta.select(band).add(delta.select('cf_t0')).subtract(delta.select('mean_cg_lr')).add(delta.select('delta'))).rename(band)

    return ndvi_cg.addBands(image.metadata('system:time_start').rename('Time')).addBands(image.select(band).subtract(ndvi_cg.select(band)).rename('error_LST'))
        .addBands(lr.select(['cluster'], ['cluster_lr'])).addBands(shortTimeFrame.cluster(clusterer).rename('cluster_t1'))
        .addBands(original.select(['cluster'], ['cluster_t0']))
        .copyProperties(lr, ['LR_Size'])



}

// Perform Global correction for 'Reflectance' mode
var globalAfter = function(original, roi, band, image_strips, t1, windowBefore, windowAfter) {

    // NIR
    var global_nir = globalCorrection(original, roi, 'nir', image_strips, t1, windowBefore, windowAfter);
    // RED
    var global_red = globalCorrection(original, roi, 'red', image_strips, t1, windowBefore, windowAfter);

    //NDVI -after

    var ndvi = global_nir.addBands(global_red).normalizedDifference(['nir', 'red']).rename('NDVI')
    return ndvi.addBands(original.select(band).subtract(ndvi.select(band)).rename('error_LST')).addBands(original.select('Time'))

}

//Perform Global correction for 'Reflectance' mode
var adaptAfter = function(originalImage, roi, k, band, image_strips, t1, windowBefore, windowAfter) {

    // NIR
    var adapt_nir = adaptiveCorrection(originalImage, roi, k, 'nir', image_strips, t1, windowBefore, windowAfter);
    // RED
    var adapt_red = adaptiveCorrection(originalImage, roi, k, 'red', image_strips, t1, windowBefore, windowAfter);

    //NDVI -after
    var ndvi = adapt_nir.addBands(adapt_red).normalizedDifference(['nir', 'red']).rename('NDVI')
    return ndvi.addBands(originalImage.select(band).subtract(ndvi.select(band)).rename('error_LST')).addBands(originalImage.select('Time'))

}

// get all the imagery to correct and correct by the chosen method
var collectionRecovered = function(date_predicted, band, k, roi, windowRange, method, ndvi, windowBefore, windowAfter) {

    // imagery start date
    var startWindowDate = ee.Date(ee.Number(date_predicted.get('system:time_start')).subtract(windowRange * (1000 * 60 * 60 * 24)));
    // imagery end date
    var endWindowDate = ee.Date(ee.Number(date_predicted.get('system:time_start')).add(windowRange * (1000 * 60 * 60 * 24)));

    var collectionLandsat = collection.filterDate(startWindowDate, endWindowDate)
        .filter(ee.Filter.notEquals('system:time_start', date_predicted.get('system:time_start')))

    var correct = collectionLandsat.map(function(img) {
        // 1- Adaptive; 0-Global; -1-LR
        if (method === 1) {
            if (ndvi) {
                return adaptiveCorrection(img, roi, k, band, windowBefore, windowAfter).copyProperties(img)
            } else {
                return adaptAfter(img, roi, k, band, windowBefore, windowAfter).copyProperties(img)
            }
        }
        if (method === 0) {
            if (ndvi) {
                return globalCorrection(img, roi, band, windowBefore, windowAfter).copyProperties(img)
            } else {
                return globalAfter(img, roi, band, windowBefore, windowAfter).copyProperties(img)
            }
        }
        if (method === -1) {
            return linearRegression(img, band, roi, windowBefore, windowAfter).select(['NDVI', 'Time']).copyProperties(img)
        }
    })

    return correct;

}

/////////// Predict /////////////
var predict = function(date_predicted, band, k, coR, modis) {

    var colDiff = coR.map(function(img) {

        var time_difference = (12 * 60 * 60 * 1000) //maximum of 12 hours
        var filter = ee.Filter.maxDifference(time_difference, 'system:time_start', img.get('Time'))
        var imgModis = ee.Image(modis.filter(filter).first()).resample('bilinear').reproject({
            crs: img.select(band).projection()
        })

        return img.addBands(img.select(band).subtract(imgModis.select(band)).rename('dif'))

    })

    // Interpolation
    var before = colDiff.first().select(['Time', 'dif'])
    var after = colDiff.sort('system:time_start', false).first().select(['Time', 'dif'])

    // Compute the ratio between the image times.
    var x1 = before.select('Time').double()
    var x2 = after.select('Time').double()
    var now = ee.Image.constant(date_predicted.get('system:time_start')).double();
    var ratio = now.subtract(x1).divide(x2.subtract(x1)) // this is zero anywhere x1 = x2
    // Compute the interpolated image.
    var interp = after.subtract(before).multiply(ratio).add(before)

    var time_difference = (12 * 60 * 60 * 1000) //maximum of 12 hours
    var filter = ee.Filter.maxDifference(time_difference, 'system:time_start', date_predicted.get('system:time_start'))

    var lst_predicted = ee.Image(modis.filter(filter).first())
        .resample('bilinear').reproject({
            crs: before.select('dif').projection()
        })
        .select(band)
        .add(interp.select('dif').rename(band))

    return lst_predicted.addBands(date_predicted.select(band).subtract(lst_predicted.select(band)).rename('error_LST'));
}

var rsquared = function(original, estimated, band, roi) {

    //sse =∑ni (yi−y^i)2
    var sse = original.select(band).subtract(estimated.select(band)).pow(2).rename(band);
    var Sum_sse = sse.reduceRegion(ee.Reducer.sum(), roi, 30, original.select(band).projection()).get(band)

    var mean = original.reduceRegion(ee.Reducer.mean(), roi, 30, original.select(band).projection()).get('NDVI')

    var ssto = original.select(band).subtract(ee.Number(mean)).pow(2);

    var Sum_ssto = ssto.reduceRegion(ee.Reducer.sum(), roi, 30, original.select(band).projection()).get(band)

    // Rsquare = 1 - (first sum/second sum)
    return ee.Number(1).subtract(ee.Number(Sum_sse).divide(ee.Number(Sum_ssto)))

}

var histograms = function(recovered, original, roi, roi_name) {

    recovered = recovered.addBands(original.select(['NDVI'], ['NDVI_LSTPredicted']))

    var options = {
        title: 'NDVI histogram in ' + roi_name,
        fontSize: 10,
        hAxis: {
            title: 'NDVI'
        },
        vAxis: {
            title: 'count of NDVI'
        },
        series: {
            0: {
                color: '#d8b365'
            },
            1: {
                color: '#5ab4ac'
            }
        }
    };

    var histogram = ui.Chart.image.histogram(recovered.select(['NDVI', 'NDVI_LSTPredicted']), roi, 30)
        .setSeriesNames(['NDVI_Predicted', 'NDVI_LST'])
        .setOptions(options);

    // Display the histogram.
    print(histogram);


}

var scatters = function(original, recovered, roi, roi_name) {


    var image = original.select('NDVI').addBands(recovered.select(['NDVI'], ['NDVI_Predicted'])).reduceRegion(ee.Reducer.toList(), roi, 30)

    var x = image.get('NDVI');
    var y = image.get('NDVI_Predicted');
    var chart = ui.Chart.array.values(y, 0, x);

    chart = chart.setSeriesNames(['NDVI'])
    chart = chart.setOptions({
        title: 'NDVI ScatterPlot in ' + roi_name,
        hAxis: {
            title: "LST Original NDVI",
            viewWindow: {
                min: 0,
                max: 1
            }
        },
        vAxis: {
            title: 'Predicted NDVI',
            viewWindow: {
                min: 0,
                max: 1
            }
        },
        pointSize: 1,
        series: {
            0: {
                color: '#d8b365'
            },
            1: {
                color: '#5ab4ac'
            }
        },
    });
    print(chart);


}


//////////////////// Visualization //////////////////////////////////

var app = {};
/** Creates the app helper functions. */
app.createHelpers = function() {
    /**
     * Enables or disables loading mode.
     * @param {boolean} enabled Whether loading mode is enabled.
     */
    app.setLoadingMode = function(enabled) {
        // Set the loading label visibility to the enabled mode.
        //app.filters.loadingLabel.style().set('shown', enabled);
        // Set each of the widgets to the given enabled mode.
        var loadDependentWidgets = [
            app.roi.select,
            app.filters.startDate,
            app.filters.endDate,
            app.method.method.select,
            app.method.K.select,
            app.filters.applyButton,
            app.filters.applyRecovered,
            app.picker.imageToCorrect.select,
            app.picker.recovered.select,
            app.picker.lr.select
        ];
        loadDependentWidgets.forEach(function(widget) {
            widget.setDisabled(enabled);
        });
    };

    /** Applies the selection filters currently selected in the UI. */
    app.applyDates = function() {
        app.setLoadingMode(true);

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi

        var filteredL8 = ee.ImageCollection(L8).filterBounds(roi)
        var filteredL7 = ee.ImageCollection(L7).filterBounds(roi)

        var start = app.filters.startDate.getValue();
        if (start) start = ee.Date(start);
        var end = app.filters.endDate.getValue();
        if (end) end = ee.Date(end);
        if (start) filteredL7 = filteredL7.filterDate(start, end).map(maskClouds).select(['B4', 'B3'], ['nir', 'red']);
        if (start) filteredL8 = filteredL8.filterDate(start, end).map(maskClouds).select(['B5', 'B4'], ['nir', 'red']);

        var filtered = filteredL7.merge(filteredL8)
            .map(addNDVI)
            .map(addTime)
            .map(function(img) {
                var clipped = img.clip(roi);
                var img_nulls = clipped.unmask(-2);
                var total_pixels = ee.Dictionary(['TOTAL_PIXELS', clipped.unmask(-2).reduceRegion(ee.Reducer.count(), roi, 30).get('NDVI')]);
                var null_pixels = ee.Dictionary(['NULL_PIXELS', img_nulls.updateMask(img_nulls.eq(-2)).reduceRegion(ee.Reducer.count(), roi, 30).get('NDVI')]);
                var quality_pixels = ee.Dictionary(['QUALITY_PIXELS', clipped.reduceRegion(ee.Reducer.count(), roi, 30).get('NDVI')]);
                var cloud_cover_region = ee.Dictionary(['CLOUD_COVER_REGION', ee.Number(null_pixels.get('NULL_PIXELS')).divide(total_pixels.get('TOTAL_PIXELS')).multiply(100)]);
                var dict = total_pixels.combine(null_pixels).combine(quality_pixels).combine(cloud_cover_region);
                var image = clipped.set(dict);

                return image.addBands(ee.Image.constant(null_pixels.get('NULL_PIXELS')).rename('Null_Pixels'))
                    .addBands(ee.Image.constant(cloud_cover_region.get('CLOUD_COVER_REGION')).rename('Cloud_Cover_Region'))
                    .addBands(ee.Image.constant(quality_pixels.get('QUALITY_PIXELS')).rename('Quality_Pixels'))
                    .addBands(ee.Image.pixelLonLat()).set('Date', ee.String(image.get('system:index')).slice(14, 22))
            })

        // Get the list of computed ids.
        var computedIdsL8 = filtered.filterMetadata('SATELLITE', 'equals', 'LANDSAT_8')
            .limit(app.IMAGE_COUNT_LIMIT)
            .reduceColumns(ee.Reducer.toList(), ['system:index'])
            .get('list');

        computedIdsL8.evaluate(function(ids) {
            // Update the image picker with the given list of ids.
            app.setLoadingMode(false);
            app.picker.imageToCorrect.select.items().reset(ids);
            app.picker.imageToCorrect.select.setValue(app.picker.imageToCorrect.select.items().get(0));
        });

        var modis = ee.ImageCollection(Modis)
            .select(['sur_refl_b02', 'sur_refl_b01'], ['nir', 'red'])
            .map(addTime)
            .filterDate(start, end)
            .filterBounds(roi)
            .map(addNDVI)
            .map(function(img) {
                return img.addBands(ee.Image.pixelLonLat())
            })

        collectionModis = modis;
        collection = filtered;
    };

    app.applyLR = function() {

        app.setLoadingMode(true);
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var time_img = imageToCorrect.get('system:time_start')
        var time = 'system:time_start'

        // imagery start date
        var startWindowDate = ee.Date(ee.Number(time_img).subtract(app.LR.daysBefore.getValue() * (1000 * 60 * 60 * 24)));
        // imagery end date
        var endWindowDate = ee.Date(ee.Number(time_img).add(app.LR.daysAfter.getValue() * (1000 * 60 * 60 * 24)));

        var colLR = collection.filterDate(startWindowDate, endWindowDate)
            .filter(ee.Filter.notEquals(time, time_img))
            .filterMetadata('CLOUD_COVER_REGION', 'less_than', 5)

        var computedIds = colLR
            .limit(app.IMAGE_COUNT_LIMIT)
            .reduceColumns(ee.Reducer.toList(), ['system:index'])
            .get('list');

        computedIds.evaluate(function(ids) {
            // Update the image picker with the given list of ids.
            app.setLoadingMode(false);
            app.picker.lr.select.items().reset(ids);
            app.picker.lr.select.setValue(app.picker.lr.select.items().get(0));
        });

    };

    app.applyRecovered = function() {

        app.setLoadingMode(true);
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var windowRange = app.filters.windowRange.getValue();
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var method = app.METHOD_OPTIONS[app.method.method.select.getValue()].value;
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var col = collectionRecovered(imageToCorrect, 'NDVI', k, roi, windowRange, method, app.method.applyNDVI.getValue(), windowBefore, windowAfter);

        var computedIds = col.limit(app.IMAGE_COUNT_LIMIT)
            .reduceColumns(ee.Reducer.toList(), ['system:index'])
            .get('list');

        computedIds.evaluate(function(ids) {
            app.setLoadingMode(false);
            app.picker.recovered.select.items().reset(ids);
            app.picker.recovered.select.setValue(app.picker.recovered.select.items().get(0));
        });

        collectionR = col;
        print(collectionR)

    }

    // Creates and styles 1 row of the legend.
    var addColor = function(color, name) {
        // Create the label that is actually the colored box.
        var colorBox = ui.Label({
            style: {
                backgroundColor: color,
                // Use padding to give the box height and width.
                padding: '9px',
                margin: '1 0 4px 0'
            }
        });

        // Create the label filled with the description text.
        var description = ui.Label({
            value: name,
            style: {
                margin: '0 0 4px 6px'
            }
        });

        return ui.Panel({
            widgets: [colorBox, description],
            layout: ui.Panel.Layout.Flow('horizontal')
        });
    };

    app.refreshMapOriginalPredicted = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))


        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())

        var NAMES = [
            'Original Image',
            'Predicted Image'
        ];

        // Create a map for each visualization option.
        var maps = [];
        NAMES.forEach(function(name, index) {
            var map = ui.Map();
            map.add(ui.Label(name));
            map.centerObject(roi, 12);
            map.addLayer(roi, {}, 'ROI');
            map.setControlVisibility(false);
            maps.push(map);
        });

        var linker = ui.Map.Linker(maps);

        // Enable zooming on the top-left map.
        maps[0].setControlVisibility({
            zoomControl: true
        });

        // Show the scale (e.g. '500m') on the bottom-right map.
        maps[1].setControlVisibility({
            scaleControl: true
        });

        // Create a grid of maps.
        var mapGrid = ui.Panel(
            [
                ui.Panel([maps[0], maps[1]], null, {
                    stretch: 'both'
                })
            ],
            ui.Panel.Layout.Flow('vertical'), {
                stretch: 'both'
            }
        );

        if (imageToCorrect) {
            // If an image id is found, create an image.
            maps[0].addLayer(imageToCorrect.select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            maps[0].addLayer(imageToCorrect.select('NDVI').unmask(-2).eq(-2).updateMask(imageToCorrect.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(imageToCorrect.select('NDVI').reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original')
            var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis)
            maps[1].addLayer(predicted.select('NDVI').clip(roi), ndviParams, 'predicted');
            maps[1].addLayer(predicted.select('NDVI').unmask(-2).eq(-2).updateMask(predicted.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');

        }

        ui.root.widgets().set(0, mapGrid);

    };

    app.refreshMapOriginalRecovered = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))


        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi;
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();

        var recovered = ee.Image(collectionR.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())
        var original = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())

        var NAMES = [
            'Original Image',
            'Recovered Image'
        ];

        // Create a map for each visualization option.
        var maps = [];
        NAMES.forEach(function(name, index) {
            var map = ui.Map();
            map.add(ui.Label(name));
            map.centerObject(roi, 12);
            map.addLayer(roi, {}, 'ROI');
            map.setControlVisibility(false);
            maps.push(map);
        });

        var linker = ui.Map.Linker(maps);

        // Enable zooming on the top-left map.
        maps[0].setControlVisibility({
            zoomControl: true
        });

        // Show the scale (e.g. '500m') on the bottom-right map.
        maps[1].setControlVisibility({
            scaleControl: true
        });

        // Create a grid of maps.
        var mapGrid = ui.Panel(
            [
                ui.Panel([maps[0], maps[1]], null, {
                    stretch: 'both'
                })
            ],
            ui.Panel.Layout.Flow('vertical'), {
                stretch: 'both'
            }
        );

        if (original && recovered) {
            maps[0].addLayer(original.select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            maps[0].addLayer(original.select('NDVI').unmask(-2).eq(-2).updateMask(original.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(original.select('NDVI').reduceRegion(ee.Reducer.mean(), roi, 30), 'Original')

            maps[1].addLayer(recovered.select('NDVI').clip(roi), ndviParams, 'Recovered');
            maps[1].addLayer(recovered.select('NDVI').unmask(-2).eq(-2).updateMask(recovered.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(recovered.select('NDVI').reduceRegion(ee.Reducer.mean(), roi, 30), 'Recovered')
            maps[1].add(legend)
        }

        ui.root.widgets().set(0, mapGrid);

    };

    app.refreshMapPredicted = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'));

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Predicted Image'))
        map.add(legend)
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())

        var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis)

        map.addLayer(predicted.select('NDVI').clip(roi), ndviParams, 'predicted');
        print(predicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean predicted');

    };

    app.refreshMapOriginal = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Original LST Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var original = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())

        map.addLayer(original.select('NDVI').clip(roi), ndviParams, 'original');
        print(original.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original LST')
    };

    app.refreshMapRecovered = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Original LST Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var recovered = ee.Image(collectionR.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())

        map.addLayer(recovered.select('NDVI').clip(roi), ndviParams, 'recovered');
        print(recovered.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean recovered LST')
    };

    app.refreshMapImagesLR = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Image from LR collection'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var original = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.lr.select.getValue()).first())

        map.addLayer(original.select('NDVI').clip(roi), ndviParams, 'lr collection image');
        print(original.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean image from LR collection')
    };

    app.refreshMapOriginalLR = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))


        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi;
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var original = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())


        var NAMES = [
            'Original Image',
            'After LR Image'
        ];

        // Create a map for each visualization option.
        var maps = [];
        NAMES.forEach(function(name, index) {
            var map = ui.Map();
            map.add(ui.Label(name));
            map.centerObject(roi, 12);
            map.addLayer(roi, {}, 'ROI');
            map.setControlVisibility(false);
            maps.push(map);
        });

        var linker = ui.Map.Linker(maps);

        // Enable zooming on the top-left map.
        maps[0].setControlVisibility({
            zoomControl: true
        });

        // Show the scale (e.g. '500m') on the bottom-right map.
        maps[1].setControlVisibility({
            scaleControl: true
        });

        // Create a grid of maps.
        var mapGrid = ui.Panel(
            [
                ui.Panel([maps[0], maps[1]], null, {
                    stretch: 'both'
                })
            ],
            ui.Panel.Layout.Flow('vertical'), {
                stretch: 'both'
            }
        );

        if (original) {
            // If an image id is found, create an image.
            maps[0].addLayer(original.select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            maps[0].addLayer(original.select('NDVI').unmask(-2).eq(-2).updateMask(original.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(original.select('NDVI').reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original')

            var lr = linearRegression(original.select('NDVI'), 'NDVI', roi, windowBefore, windowAfter);
            maps[1].addLayer(lr.select('NDVI').clip(roi), ndviParams, 'LRImage');
            maps[1].addLayer(lr.select('NDVI').unmask(-2).eq(-2).updateMask(lr.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            maps[1].add(legend)
            print(lr.select('NDVI').reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean LR')
        }

        ui.root.widgets().set(0, mapGrid);

    };

    app.refreshMapClusters = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#ffffe5', '0.2')).add(addColor('#f7fcb9', '')).add(addColor('#d9f0a3', '')).add(addColor('#addd8e', ''))
            .add(addColor('#78c679', '')).add(addColor('#41ab5d', '')).add(addColor('#238443', '')).add(addColor('#006837', '')).add(addColor('#004529', '0.8'))


        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi;
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var recovered = ee.Image(collectionR.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())

        var NAMES = [
            'Clusters ImageStripped T0',
            'Clusters ImageT1',
            'Clusters Image LR'
        ];

        // Create a map for each visualization option.
        var maps = [];
        NAMES.forEach(function(name, index) {
            var map = ui.Map();
            map.add(ui.Label(name));
            map.centerObject(roi, 12);
            map.addLayer(roi, {}, 'ROI');
            map.setControlVisibility(false);
            maps.push(map);
        });

        var linker = ui.Map.Linker(maps);

        // Enable zooming on the top-left map.
        maps[0].setControlVisibility({
            zoomControl: true
        });

        // Show the scale (e.g. '500m') on the bottom-right map.
        maps[1].setControlVisibility({
            scaleControl: true
        });

        // Create a grid of maps.
        var mapGrid = ui.Panel(
            [
                ui.Panel([maps[0], maps[1], maps[2]], null, {
                    stretch: 'both'
                })
            ],
            ui.Panel.Layout.Flow('vertical'), {
                stretch: 'both'
            }
        );

        if (recovered.select('cluster_t0')) {

            maps[0].addLayer(recovered.select('cluster_t0').randomVisualizer().clip(roi), {}, 'clusters T0');
            maps[1].addLayer(recovered.select('cluster_t1').randomVisualizer().clip(roi), {}, 'clusters T1');
            maps[2].addLayer(recovered.select('cluster_lr').randomVisualizer().clip(roi), {}, 'clusters LR');
        } else {
            print('Impossible to Show Clusters. Select method Adaptive or NDVI from the beginning')
        }

        ui.root.widgets().set(0, mapGrid);
    };

    app.refreshMapError = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#fc8d59', '-0.2')).add(addColor('#ffffbf', '')).add(addColor('#91cf60', '0.2'))

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Error LST'))
        map.add(legend)
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())

        var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis)

        map.addLayer(predicted.select('error_LST').clip(roi), errorParams, 'error');
        print(predicted.select('error_LST').reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean error LST')
    };




    /// STATS ///////
    app.showHistograms = function() {

        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var k_name = app.K_OPTIONS[app.method.K.select.getValue()].description;


        var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis)

        var image = predicted.select('NDVI').mask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
        imageToCorrect = imageToCorrect.select('NDVI').mask(predicted.select('NDVI').unmask(-2).neq(-2))
        histograms(image, imageToCorrect, roi, roi_name);

        //})

        /*   var adapt4 = adaptiveCorrection(imageToCorrect,roi,4,'NDVI',stripped, t1, windowBefore, windowAfter).select(['NDVI'],['NDVI_K4'])
    var adapt6 = adaptiveCorrection(imageToCorrect,roi,6,'NDVI',stripped, t1, windowBefore, windowAfter).select(['NDVI'],['NDVI_K6'])
    var adapt8 = adaptiveCorrection(imageToCorrect,roi,8,'NDVI',stripped, t1, windowBefore, windowAfter).select(['NDVI'],['NDVI_K8'])
    var image = adapt4.addBands(adapt6.select('NDVI_K6')).addBands(adapt8.select('NDVI_K8')).mask(stripped.select('NDVI').unmask(-2).eq(-2))
    print(image)

    var histogramW = ui.Chart.image.histogram(image.select(['NDVI_K4', 'NDVI_K6','NDVI_K8']), roi, 30)
    .setSeriesNames(['NDVI_K4', 'NDVI_K6','NDVI_K8'])
    .setOptions({
    title: 'NDVI histogram in ' + roi_name +' Using Adaptive Method',
    fontSize: 10,
    hAxis: {title: 'NDVI'},
    vAxis: {title: 'count of NDVI'}
    });
print(histogramW);
*/


    }

    app.showScatters = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var k_name = app.K_OPTIONS[app.method.K.select.getValue()].description;

        var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis)

        var image = predicted.select('NDVI').mask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
        imageToCorrect = imageToCorrect.select('NDVI').mask(predicted.select('NDVI').unmask(-2).neq(-2))

        scatters(imageToCorrect, image, roi, roi_name)


    }

    app.metrics = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var windowBefore = app.LR.daysBefore.getValue();
        var windowAfter = app.LR.daysAfter.getValue();
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
        var k_name = app.K_OPTIONS[app.method.K.select.getValue()].description;

        var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis);

        var image = predicted.select('NDVI').mask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
        imageToCorrect = imageToCorrect.select('NDVI').mask(predicted.select('NDVI').unmask(-2).neq(-2))

        print(rsquared(imageToCorrect, image, 'NDVI', roi), 'R squared')

        var error = imageToCorrect.select('NDVI').subtract(image.select('NDVI')).rename('error')
        var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
        var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
        var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])

        var dict = mae.combine(rmse).combine(error_dict);
        print(dict, 'Stats')
        print(image.addBands(imageToCorrect.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30), 'pearsons Correlation');

    }

    app.regionStats = function() {

        var ROI = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var ROI_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR

        var stats = Chart.image.series(collection.select(['Null_Pixels', 'Quality_Pixels']), ROI, ee.Reducer.first(), 30)
            .setOptions({
                title: 'Nº Pixels in ' + ROI_name,
                vAxis: {
                    title: 'Nº Pixels'
                },
                hAxis: {
                    title: 'Date'
                },
            });

        print(stats)

        var stats_percent = Chart.image.series(collection.select(['Cloud_Cover_Region']), ROI, ee.Reducer.first(), 30)
            .setOptions({
                title: '% Cloud Cover in ' + ROI_name,
                vAxis: {
                    title: '% Cloud Pixels'
                },
                hAxis: {
                    title: 'Date'
                },
            });

        print(stats_percent)

        var stats_ndvi = Chart.image.series(collection.select(['NDVI']), ROI, ee.Reducer.mean(), 30)
            .setOptions({
                title: 'NDVI Mean in ' + ROI_name + ' in All Dates',
                vAxis: {
                    title: 'NDVI Mean'
                },
                hAxis: {
                    title: 'Date'
                },
            });

        print(stats_ndvi)

        var stats_ndvi_2 = Chart.image.series(collection.filterMetadata('CLOUD_COVER_REGION', 'less_than', 1).select(['NDVI']), ROI, ee.Reducer.mean(), 30)
            .setOptions({
                title: 'NDVI Mean in ' + ROI_name + ' with Cloud_Cover_Region < 1%',
                vAxis: {
                    title: 'NDVI Mean'
                },
                hAxis: {
                    title: 'Date'
                },
            });

        print(stats_ndvi_2)

    }

}

////////// Export ////////////////////

app.exportOriginal = function() {

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
    var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
    var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())


    // Export Original Image
    Export.image.toDrive({
        image: imageToCorrect.select('NDVI').clip(roi).visualize(ndviParams),
        description: 'LST_OriginalImage-' + app.picker.imageToCorrect.select.getValue(),
        scale: 30,
        region: roi
    });

}

app.exportPredicted = function() {

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
    var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
    var windowBefore = app.LR.daysBefore.getValue();
    var windowAfter = app.LR.daysAfter.getValue();
    var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
    var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis);


    Export.image.toDrive({
        image: predicted.select('NDVI').clip(roi).visualize(ndviParams),
        description: 'Predicted-' + app.picker.imageToCorrect.select.getValue(),
        scale: 30,
        region: roi
    });

}

app.exportError = function() {

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
    var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
    var windowBefore = app.LR.daysBefore.getValue();
    var windowAfter = app.LR.daysAfter.getValue();
    var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
    var predicted = predict(imageToCorrect, 'NDVI', k, collectionR, collectionModis);

    Export.image.toDrive({
        image: predicted.select('error_LST').clip(roi).visualize(errorParams),
        description: 'ErrorPredicted-' + app.picker.imageToCorrect.select.getValue(),
        scale: 30,
        region: roi
    });

}

app.exportRecovered = function() {

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
    var k = app.K_OPTIONS[app.method.K.select.getValue()].k;
    var windowBefore = app.LR.daysBefore.getValue();
    var windowAfter = app.LR.daysAfter.getValue();
    var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
    var recovered = ee.Image(collectionR.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())

    Export.image.toDrive({
        image: recovered.select('NDVI').clip(roi).visualize(ndviParams),
        description: 'recovered-' + app.picker.recovered.select.getValue(),
        scale: 30,
        region: roi
    });

}

app.exportTables = function() {

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
    var method = app.METHOD_OPTIONS[app.method.method.select.getValue()].value;
    var method_name = app.METHOD_OPTIONS[app.method.method.select.getValue()].description;
    var signal = app.method.applyNDVI.getValue();
    var windowRange = app.filters.windowRange.getValue();
    var windowBefore = app.LR.daysBefore.getValue();
    var windowAfter = app.LR.daysAfter.getValue();
    var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
    var date = imageToCorrect.get('Date').getInfo()
    var k_method = app.K_OPTIONS[app.method.K.select.getValue()].k;

    var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
    var roi_abre = app.ROI_OPTIONS[app.roi.select.getValue()].abre

    var k = [4, 6, 8];
    var dici = ee.Dictionary();
    var features = ee.List([]);
    var featureswatermasked = ee.List([]);
    var featuresroadsmasked = ee.List([]);
    var featuresbuildingsmasked = ee.List([]);
    var featureswrmasked = ee.List([]);
    var featureswbmasked = ee.List([]);
    var featuresrbmasked = ee.List([]);
    var featuresallmasked = ee.List([]);
    
    var calc_Metrics_adaptive = function(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col) {
        var r = ee.Number(rsquared(imageToCorrect, image, 'NDVI', roi))

        // Original - Predicted
        var r_dict = ee.Dictionary(['R^2', r])
        var method = ee.Dictionary(['method', 'adaptive'])
        var signal = ee.Dictionary(['signal', 1])
        var k_dict = ee.Dictionary(['K', nameK])
        var LRA_dict = ee.Dictionary(['LRA', windowAfter])
        var LRB_dict = ee.Dictionary(['LRB', windowBefore])
        var WR_dict = ee.Dictionary(['WR', windowRange])
        var pairs = ee.Dictionary(['Pairs', col.size()])
        var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
        var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
        var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
        var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
        var corr = image.addBands(imageToCorrect.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
        var dict = method.combine(signal).combine(k_dict).combine(LRA_dict).combine(LRB_dict).combine(WR_dict).combine(pairs)
            .combine(r_dict).combine(mae).combine(rmse).combine(error_dict).combine(corr)
      
        return ee.Feature(null, dict);
    }
    
    var calc_Metrics_global = function(date_predicted, image, windowBefore, windowAfter, windowRange, col) {
        var r = ee.Number(rsquared(imageToCorrect, image, 'NDVI', roi))

        // Original - Predicted
        var r_dict = ee.Dictionary(['R^2', r])
        var method = ee.Dictionary(['method', 'global'])
        var signal = ee.Dictionary(['signal', 0])
        var LRA_dict = ee.Dictionary(['LRA', windowAfter])
        var LRB_dict = ee.Dictionary(['LRB', windowBefore])
        var WR_dict = ee.Dictionary(['WR', windowRange])
        var pairs = ee.Dictionary(['Pairs', col.size()])
        var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
        var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
        var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
        var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
        var corr = image.addBands(imageToCorrect.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
        var dict = method.combine(signal).combine(LRA_dict).combine(LRB_dict).combine(WR_dict).combine(pairs)
            .combine(r_dict).combine(mae).combine(rmse).combine(error_dict).combine(corr)
      
        return ee.Feature(null, dict);
    }
    
    var makePredictionNDVIGlobal = function(windowBefore, windowAfter, windowRange, col, imageToCorrect) {
        // NDVI
        var imageNDVI = predict(imageToCorrect, 'NDVI', k_method, col, collectionModis);
        
        var image = imageNDVI.select('NDVI').clip(roi).mask(imageToCorrect.select('NDVI').clip(roi).unmask(-2).neq(-2))
        var date_predicted = imageToCorrect.select('NDVI').clip(roi).mask(imageToCorrect.select('NDVI').clip(roi).unmask(-2).neq(-2))

        var ft = calc_Metrics_global(date_predicted, image, windowBefore, windowAfter, windowRange, col)
        return ee.List([ft]);
    }
    
    var createBuffer = function(feature){ return feature.buffer(20); }
    
    if(method === 1) {
      if(app.method.applyNDVI.getValue())
      {
        // Adaptive - NDVI
        k.forEach(function(nameK, indexK) {
            var col = collectionRecovered(imageToCorrect, 'NDVI', k[indexK], roi, windowRange, 1, true, windowBefore, windowAfter);
            var imageNDVI = predict(imageToCorrect, 'NDVI', k[indexK], col, collectionModis);
            var NDVIwaterMask = imageNDVI.clip(osmwater).unmask().not();
            var NDVIroadMask = imageNDVI.clip(osmroads.map(createBuffer)).unmask().not();
            var NDVIbuildingMask = imageNDVI.clip(osmbuildings).unmask().not();
            

            // NO MASK
            var image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
            var date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
            var ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            features = features.cat(ee.List([ft]))
            
            // WATER MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswatermasked = featureswatermasked.cat(ee.List([ft]))

            // ROAD MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresroadsmasked = featuresroadsmasked.cat(ee.List([ft]))

            // BUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresbuildingsmasked = featuresbuildingsmasked.cat(ee.List([ft]))

            // WATERROAD MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswrmasked = featureswrmasked.cat(ee.List([ft]))
            
            // WATERBUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswbmasked = featureswbmasked.cat(ee.List([ft]))
            
            // ROADBUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresrbmasked = featuresrbmasked.cat(ee.List([ft]))
            
            // ALL MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresallmasked = featuresallmasked.cat(ee.List([ft]))
        })
  
      } else {
    
        // Adaptive - Reflectance
        k.forEach(function(nameK, indexK) {
            var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
            var col = collectionRecovered(imageToCorrect, 'NDVI', k[indexK], roi, windowRange, 1, true, windowBefore, windowAfter);
            var imageNDVI = predict(imageToCorrect, 'NDVI', k[indexK], col, collectionModis);
            var NDVIwaterMask = imageNDVI.clip(osmwater).unmask().not();
            var NDVIroadMask = imageNDVI.clip(osmroads.map(createBuffer)).unmask().not();
            var NDVIbuildingMask = imageNDVI.clip(osmbuildings).unmask().not();
            

            // NO MASK
            var image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
            var date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2))
            var ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            features = features.cat(ee.List([ft]))
            
            // WATER MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswatermasked = featureswatermasked.cat(ee.List([ft]))

            // ROAD MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresroadsmasked = featuresroadsmasked.cat(ee.List([ft]))

            // BUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresbuildingsmasked = featuresbuildingsmasked.cat(ee.List([ft]))

            // WATERROAD MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswrmasked = featureswrmasked.cat(ee.List([ft]))
            
            // WATERBUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featureswbmasked = featureswbmasked.cat(ee.List([ft]))
            
            // ROADBUILDING MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresrbmasked = featuresrbmasked.cat(ee.List([ft]))
            
            // ALL MASK
            image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            date_predicted = imageToCorrect.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(imageToCorrect.select('NDVI').unmask(-2).neq(-2)).updateMask(NDVIwaterMask).updateMask(NDVIroadMask).updateMask(NDVIbuildingMask)
            ft = calc_Metrics_adaptive(date_predicted, image, windowBefore, windowAfter, windowRange, nameK, col)
            featuresallmasked = featuresallmasked.cat(ee.List([ft]))
    
        })
      }
    }
    else if(method === 0)
    {
      if(app.method.applyNDVI.getValue())
      {
        // Global - NDVI
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var col = collectionRecovered(imageToCorrect, 'NDVI', k_method, roi, windowRange, 0, true, windowBefore, windowAfter);
        var NDVIwaterMask = imageToCorrect.clip(osmwater).unmask().not();
        var NDVIroadMask = imageToCorrect.clip(osmroads.map(createBuffer)).unmask().not();
        var NDVIbuildingMask = imageToCorrect.clip(osmbuildings).unmask().not();

        // NO MASK
        var ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect)
        features = features.cat(ft)

        // WATER MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIwaterMask.select('NDVI')))
        featureswatermasked = featureswatermasked.cat(ft)

        // ROAD MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIroadMask.select('NDVI')))
        featuresroadsmasked = featuresroadsmasked.cat(ft)

        // BUILDING MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIbuildingMask.select('NDVI')))
        featuresbuildingsmasked = featuresbuildingsmasked.cat(ft)

        // WATERROAD MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIwaterMask.select('NDVI')).mask(NDVIroadMask.select('NDVI')))
        featureswrmasked = featureswrmasked.cat(ft)
        
        // WATERBUILDING MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIwaterMask.select('NDVI')).mask(NDVIbuildingMask.select('NDVI')))
        featureswbmasked = featureswbmasked.cat(ft)
        
        // ROADBUILDING MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIroadMask.select('NDVI')).mask(NDVIbuildingMask.select('NDVI')))
        featuresrbmasked = featuresrbmasked.cat(ft)
        
        // ALL MASK
        ft = makePredictionNDVIGlobal(windowBefore, windowAfter, windowRange, col, imageToCorrect.mask(NDVIwaterMask.select('NDVI')).mask(NDVIroadMask.select('NDVI')).mask(NDVIbuildingMask.select('NDVI')))
        featuresallmasked = featuresallmasked.cat(ft)

      } else {
        // Global - Reflectance
        var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
        var col = collectionRecovered(imageToCorrect, 'NDVI', k_method, roi, windowRange, 0, false, windowBefore, windowAfter);
        var image = predict(imageToCorrect, 'NDVI', k_method, col, collectionModis);
        image = image.select('NDVI').clip(roi).mask(imageToCorrect.select('NDVI').clip(roi).unmask(-2).neq(-2))
        var date_predicted = imageToCorrect.select('NDVI').clip(roi).mask(image.select('NDVI').clip(roi).unmask(-2).neq(-2))
    
        var r = ee.Number(rsquared(imageToCorrect, image, 'NDVI', roi))
    
        // Original - Predicted
        var r_dict = ee.Dictionary(['R^2', r])
        var method = ee.Dictionary(['method', 'global'])
        var signal = ee.Dictionary(['signal', 0])
        var LRA_dict = ee.Dictionary(['LRA', windowAfter])
        var LRB_dict = ee.Dictionary(['LRB', windowBefore])
        var WR_dict = ee.Dictionary(['WR', windowRange])
        var pairs = ee.Dictionary(['Pairs', col.size()])
        var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
        var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
        var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
        var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
        var corr = image.addBands(imageToCorrect.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
        var dict = method.combine(signal).combine(LRA_dict).combine(LRB_dict).combine(WR_dict).combine(pairs)
            .combine(r_dict).combine(mae).combine(rmse).combine(error_dict).combine(corr)
    
        var ft = ee.Feature(null, dict)
        features = features.cat(ee.List([ft]))
      }
    }
    else{
      // LR
      var imageToCorrect = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
      var col = collectionRecovered(imageToCorrect, 'NDVI', k_method, roi, windowRange, -1, true, windowBefore, windowAfter);
      var image = predict(imageToCorrect, 'NDVI', k_method, col, collectionModis);
      image = image.select('NDVI').clip(roi).mask(imageToCorrect.select('NDVI').clip(roi).unmask(-2).neq(-2))
      var date_predicted = imageToCorrect.select('NDVI').clip(roi).mask(image.select('NDVI').clip(roi).unmask(-2).neq(-2))
  
      var r = ee.Number(rsquared(imageToCorrect, image, 'NDVI', roi))
  
      // Original - Predicted
      var r_dict = ee.Dictionary(['R^2', r])
      var method = ee.Dictionary(['method', 'LR'])
      var signal = ee.Dictionary(['signal', 1])
      var LRA_dict = ee.Dictionary(['LRA', windowAfter])
      var LRB_dict = ee.Dictionary(['LRB', windowBefore])
      var WR_dict = ee.Dictionary(['WR', windowRange])
      var pairs = ee.Dictionary(['Pairs', col.size()])
      var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
      var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
      var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
      var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
      var corr = image.addBands(imageToCorrect.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
      var dict = method.combine(signal).combine(LRA_dict).combine(LRB_dict).combine(WR_dict).combine(pairs)
          .combine(r_dict).combine(mae).combine(rmse).combine(error_dict).combine(corr)
  
      var ft = ee.Feature(null, dict)
      features = features.cat(ee.List([ft]))
    }
    var metrics = ee.FeatureCollection(features)

    print(metrics)

    Export.table.toDrive({
        collection: metrics,
        description: 'stair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
    
    var metrics = ee.FeatureCollection(featureswatermasked)
    Export.table.toDrive({
        collection: metrics,
        description: 'wstair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
    
    var metrics = ee.FeatureCollection(featuresroadsmasked)
    Export.table.toDrive({
        collection: metrics,
        description: 'rstair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
    
    var metrics = ee.FeatureCollection(featuresbuildingsmasked)
    Export.table.toDrive({
        collection: metrics,
        description: 'bstair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
    
    var metrics = ee.FeatureCollection(featureswrmasked)
    Export.table.toDrive({
        collection: metrics,
        description: '1stair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
            
    var metrics = ee.FeatureCollection(featureswbmasked)
    Export.table.toDrive({
        collection: metrics,
        description: '2stair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
    
    var metrics = ee.FeatureCollection(featuresrbmasked)
    Export.table.toDrive({
        collection: metrics,
        description: '3stair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
            
    var metrics = ee.FeatureCollection(featuresallmasked)
    Export.table.toDrive({
        collection: metrics,
        description: '4stair_' + roi_abre + '_R' + windowRange + '_B' + windowBefore + 'A' + windowAfter + '_' + date,
        fileFormat: 'CSV',
        folder: 'metrics2',
        selectors: ['Pairs', 'WR', 'signal', 'K', 'method', 'LRB', 'LRA', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
    })
}


app.createPanels = function() {

    app.intro = {
        panel: ui.Panel([
            ui.Label({
                value: 'STAIR Predict - Landsat/Modis',
                style: {
                    fontWeight: 'bold',
                    fontSize: '24px',
                    margin: '10px 5px',
                    shown: true
                }
            })
        ])
    };

    app.roi = {
        label: ui.Label(),
        // Create a select with a function that reacts to the "change" event.
        select: ui.Select({
            items: Object.keys(app.ROI_OPTIONS),
            onChange: function() {
                // Update the label's value with the select's description.
                var option = app.ROI_OPTIONS[app.roi.select.getValue()];
                app.roi.label.setValue(option.description);
            }
        })
    }

    app.roi.panel = ui.Panel({
        widgets: [
            ui.Label('1) Select Region of Interest', {
                fontWeight: 'bold'
            }),
            app.roi.select
        ],
        style: app.SECTION_STYLE
    })

    app.roi.select.setValue(app.roi.select.items().get(0));

    app.method = {
        applyNDVI: ui.Checkbox({
            label: 'Apply NDVI at first',
            value: true
        }),
        method: {
            label: ui.Label(),
            // Create a select with a function that reacts to the "change" event.
            select: ui.Select({
                items: Object.keys(app.METHOD_OPTIONS),
                onChange: function() {
                    // Update the label's value with the select's description.
                    var option = app.METHOD_OPTIONS[app.method.method.select.getValue()];
                    app.method.method.label.setValue(option.description);
                }
            })
        },
        K: {
            label: ui.Label(),
            // Create a select with a function that reacts to the "change" event.
            select: ui.Select({
                items: Object.keys(app.K_OPTIONS),
                onChange: function() {
                    // Update the label's value with the select's description.
                    var option = app.K_OPTIONS[app.method.K.select.getValue()];
                    app.method.K.label.setValue(option.description);
                }
            })
        }
    }

    app.method.panel = ui.Panel({
        widgets: [
            ui.Label('2) Select method to recover images', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Method', app.HELPER_TEXT_STYLE), app.method.method.select,
                ui.Label('K', app.HELPER_TEXT_STYLE), app.method.K.select,
            ], ui.Panel.Layout.flow('horizontal')), app.method.applyNDVI
        ],
        style: app.SECTION_STYLE
    })

    app.method.method.select.setValue(app.method.method.select.items().get(0));
    app.method.K.select.setValue(app.method.K.select.items().get(0));

    app.filters = {
        startDate: ui.Textbox('YYYY-MM-DD', '2017-01-01'),
        endDate: ui.Textbox('YYYY-MM-DD', '2017-12-31'),
        applyButton: ui.Button('Apply filters', app.applyDates),
        applyRecovered: ui.Button('Get recovered images', app.applyRecovered),
        loadingLabel: ui.Label({
            value: 'Loading...',
            style: {
                stretch: 'vertical',
                color: 'gray',
                shown: false
            }
        }),
        windowRange: ui.Textbox({
            placeholder: 'NN',
            value: '10',
            style: {
                maxWidth: '50px'
            }
        })
    };

    app.filters.panel = ui.Panel({
        widgets: [
            ui.Label('3) Select Dates of Interest', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Panel([
                    ui.Label('Start date', app.HELPER_TEXT_STYLE), app.filters.startDate
                ], ui.Panel.Layout.flow('vertical')),
                ui.Panel([
                    ui.Label('End date', app.HELPER_TEXT_STYLE), app.filters.endDate
                ], ui.Panel.Layout.flow('vertical'))
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                app.filters.applyButton,
                app.filters.loadingLabel
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

    app.picker = {
        // Create a select with a function that reacts to the "change" event.
        imageToCorrect: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.imageToCorrect.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.imageToCorrect.label.setValue('Cloud Cover Region: ' + result.toFixed(2))
                    })
                }
            })

        },
        recovered: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.recovered.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.recovered.label.setValue('Cloud Cover Region: ' + result.toFixed(2))
                    })
                }
            })

        },
        lr: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.lr.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.lr.label.setValue('Cloud Cover Region: ' + result.toFixed(2))
                    })
                }
            })
        }
    }

    app.LR = {
        daysBefore: ui.Textbox({
            placeholder: 'NN',
            value: '10',
            style: {
                maxWidth: '50px'
            }
        }),
        daysAfter: ui.Textbox({
            placeholder: 'NN',
            value: '10',
            style: {
                maxWidth: '50px'
            }
        })
    };

    app.picker.panel = ui.Panel({
        widgets: [
            ui.Label('4) Inputs to Predict', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Image To Predict', app.HELPER_TEXT_STYLE),
                ui.Panel([
                    app.picker.imageToCorrect.select,
                    app.picker.imageToCorrect.label
                ], ui.Panel.Layout.flow('horizontal'))
            ], ui.Panel.Layout.flow('vertical')),
            ui.Panel([
                ui.Label('Window Range to Predict', app.HELPER_TEXT_STYLE), app.filters.windowRange,
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Label('Linear Regression', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Days Before', app.HELPER_TEXT_STYLE), app.LR.daysBefore,
                ui.Label('Days After', app.HELPER_TEXT_STYLE), app.LR.daysAfter,
                ui.Button('Get LR images', app.applyLR)
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Images', app.HELPER_TEXT_STYLE), app.picker.lr.select,
                app.picker.lr.label
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Label('5) Recovered Images', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                app.filters.applyRecovered
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Recovered Images', app.HELPER_TEXT_STYLE), app.picker.recovered.select,
                app.picker.recovered.label,
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });



    app.LR.panel = ui.Panel({
        widgets: [
            ui.Label('Linear Regression', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Days Before', app.HELPER_TEXT_STYLE), app.LR.daysBefore,
                ui.Label('Days After', app.HELPER_TEXT_STYLE), app.LR.daysAfter,
                ui.Button('Get LR images', app.applyLR)
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Images', app.HELPER_TEXT_STYLE), app.picker.lr.select,
                app.picker.lr.label
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    })

    app.show = {
        buttonOriginalPredicted: ui.Button({
            label: 'Show Original - Predicted',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalPredicted()
            }
        }),
        buttonOriginalRecovered: ui.Button({
            label: 'Show Original - Recovered',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalRecovered()
            }
        }),
        buttonOriginal: ui.Button({
            label: 'Show Original',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginal()
            }
        }),
        buttonPredicted: ui.Button({
            label: 'Show Predicted',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapPredicted()
            }
        }),
        buttonRecovered: ui.Button({
            label: 'Show Recovered',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapRecovered()
            }
        }),
        buttonClusters: ui.Button({
            label: 'Clusters Recovered',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapClusters()
            }
        }),
        buttonOriginalLR: ui.Button({
            label: 'Show Original-LR',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalLR()
            }
        }),
        buttonImagesLR: ui.Button({
            label: 'Show Images LR',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapImagesLR()
            }
        }),
        buttonError: ui.Button({
            label: 'Show Errors',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapError()
            }
        }),
        buttonHistograms: ui.Button({
            label: 'Histograms',
            // React to the button's click event.
            onClick: function() {
                app.showHistograms()
            }
        }),
        buttonScatters: ui.Button({
            label: 'Scatter Plots',
            // React to the button's click event.
            onClick: function() {
                app.showScatters()
            }
        }),
        buttonStatsRegion: ui.Button({
            label: 'Region Stats',
            // React to the button's click event.
            onClick: function() {
                app.regionStats()
            }
        }),
        buttonMetrics: ui.Button({
            label: 'Metrics',
            // React to the button's click event.
            onClick: function() {
                app.metrics()
            }
        })
    }

    app.show.panel = ui.Panel({
        widgets: [
            ui.Label('5) Visualization', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Using Method', app.HELPER_TEXT_STYLE), app.show.buttonOriginalRecovered,
                app.show.buttonRecovered
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                app.show.buttonClusters, app.show.buttonOriginalLR, app.show.buttonImagesLR
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Predict', app.HELPER_TEXT_STYLE), app.show.buttonOriginal, app.show.buttonOriginalPredicted
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                app.show.buttonPredicted, app.show.buttonError
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Others', app.HELPER_TEXT_STYLE),
                app.show.buttonMetrics, app.show.buttonStatsRegion, app.show.buttonScatters, app.show.buttonHistograms
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

    app.export = {

        buttonRecoveredExport: ui.Button({
            label: 'Export Predicted',
            // React to the button's click event.
            onClick: function() {
                app.exportPredicted();
            }
        }),
        buttonOriginalExport: ui.Button({
            label: 'Export Original',
            // React to the button's click event.
            onClick: function() {
                app.exportOriginal();
            }
        }),
        buttonErrorExport: ui.Button({
            label: 'Export Error',
            // React to the button's click event.
            onClick: function() {
                app.exportError();
            }
        }),
        buttonImagesExport: ui.Button({
            label: 'Export Recovered Images',
            // React to the button's click event.
            onClick: function() {
                app.exportRecovered();
            }
        }),
        buttonExportTable: ui.Button({
            label: 'Export Tables',
            // React to the button's click event.
            onClick: function() {
                app.exportTables();
            }
        })
    };

    app.export.panel = ui.Panel({
        widgets: [
            ui.Label('5) Start an export', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                app.export.buttonRecoveredExport,
                app.export.buttonOriginalExport,
                app.export.buttonErrorExport,
                app.export.buttonImagesExport,
                app.export.buttonExportTable
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

}

app.createConstants = function() {
    app.COLLECTION_ID = 'LANDSAT/LC08/C01/T1_SR';
    app.COLLECTION_ID_L7 = 'LANDSAT/LE07/C01/T1_SR';
    app.SECTION_STYLE = {
        margin: '20px 0 0 0'
    };
    app.HELPER_TEXT_STYLE = {
        margin: '8px 0 -3px 8px',
        fontSize: '12px',
        color: 'gray'
    };
    app.IMAGE_COUNT_LIMIT = 20;
    app.NDVI_OPTIONS = {
        'Reflectance': {
            description: 'Applying NDVI after the prediction of NIR and RED'
        },
        'NDVI': {
            description: 'Predict NDVI from te beginning'
        }
    };
    app.ROI_OPTIONS = {
        'São Cristóvão de Lafões': {
            nameR: 'Santa Cruz da Trapa e São Cristovão de Lafões',
            abre: 'SCTSCL',
            roi: saocristovao
        },
        'Benfica do Ribatejo ': {
            nameR: 'Benfira do Ribatejo',
            abre: 'BenRib',
            roi: benfica
        },
        'Santar': {
            nameR: 'Santar',
            abre: 'Santar',
            roi: santar
        },
        'Beringel': {
            nameR: 'Beringel',
            abre: 'BGL',
            roi: beringel
        }
    };
    app.K_OPTIONS = {
        '4': {
            description: '4',
            k: 4
        },
        '6': {
            description: '6',
            k: 6
        },
        '8': {
            description: '8',
            k: 8
        }
    }
    app.METHOD_OPTIONS = {
        'Adaptive': {
            description: 'Adaptive',
            value: 1
        },
        'Global': {
            description: 'Global',
            value: 0
        },
        'Linear Regression': {
            description: 'LR',
            value: -1
        }
    }
};

app.boot = function() {
    app.createConstants();
    app.createHelpers();
    app.createPanels();
    var main = ui.Panel({
        widgets: [
            app.intro.panel,
            app.roi.panel,
            app.method.panel,
            app.filters.panel,
            app.picker.panel,
            app.show.panel,
            app.export.panel
        ],
        style: {
            width: '430px',
            padding: '8px'
        }
    });
    Map.setCenter(-7.799, 39.20, 6);
    ui.root.insert(1, main);
};

app.boot();