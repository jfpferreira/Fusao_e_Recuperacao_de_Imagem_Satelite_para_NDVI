var L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    modis = ee.ImageCollection("MODIS/006/MOD09GQ"),
    osmwater = ee.FeatureCollection("users/jfpferreira/OSM_Water"),
    osmroads = ee.FeatureCollection("users/jfpferreira/OSM_Roads"),
    osmbuildings = ee.FeatureCollection("users/jfpferreira/OSM_Buildings"),
//ROI
var saocristovao = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '181621');
var santar = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '180911');
var benfica = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '140302');
var beringel = ee.FeatureCollection('users/danielahenriques16/Cont_AAD_CAOP2017').filterMetadata('Dicofre', 'equals', '020503');

// Palette
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
var collection = ee.ImageCollection(ee.Image.constant(0))

// A function to mask out cloudy pixels.
var maskCloudsModis = function(image) {
    var QA = image.select('QC_250m')
    // Make a mask to get bit 10, the internal_cloud_algorithm_flag bit.
    var bitMask = 1 << 1;
    // Return an image masking out cloudy areas.
    return image.updateMask(QA.bitwiseAnd(bitMask).eq(0))
}

var maskCloudsLandsat = function(image) {
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

var addNDVI = function(img) {
    return img.addBands(img.normalizedDifference(['nir', 'red']).rename('NDVI'))
};

var addTime = function(img) {
    return img
        // System Time
        .addBands(img.metadata('system:time_start').rename('Time'))
        // Date - Day is the Julian date
        .set('Date', ee.Date(img.get('system:time_start')).format('YYYY-MM-DD'))
        // Longitude and Latitude
        .addBands(ee.Image.pixelLonLat())
};

/*
date_pair - Landsat Image at date 1 (pair)
date_predicted - original image Landsat, used to compare the predicted image to the original
Window_pixels - Window, can be 3,5,7,9
K - use for clustering. Can be 4,6,8
band - NIR, RED or NDVI
ROI - From the 4 available
*/
var predict = function(window_pixels, k, band, roi, date_pair, date_predicted) {

    var modis_pair = ee.Image(modis
        .filterBounds(roi)
        .map(addTime)
        .filterMetadata('Date', 'equals', date_pair.get('Date'))
        .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
        .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
        crs: date_pair.select(band).projection()
    })

    var modis_predicted = ee.Image(modis
        .filterBounds(roi)
        .map(addTime)
        .filterMetadata('Date', 'equals', date_predicted.get('Date'))
        .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
        .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
        crs: date_predicted.select(band).projection()
    })

    var spectralDifference = function(img_landsat, img_modis, band) {

        return img_landsat.select(band).subtract(img_modis.select(band)).abs().rename(['SpectralDiff'])

    }

    var temporalDifference = function(img_pair, img_predicted) {

        var diff = img_predicted.select('Time').subtract(img_pair.select('Time')).abs()
        return ee.Image(diff).rename(['TempDiff']);
    }

    // Landsat_Pair with TempDiff/ SpectralDiff/ NDVI_ModisPair/ NDVI_ModisPredicted
    var lst_img = date_pair.select([band, 'latitude', 'longitude'], [band + '_LSTPair', 'latitude', 'longitude'])
        .addBands(temporalDifference(date_pair, date_predicted), ['TempDiff'])
        .addBands(spectralDifference(date_pair, modis_pair, band), ['SpectralDiff'])
        .addBands(modis_pair.select([band], [band + '_ModisPair']))
        .addBands(modis_predicted.select([band], [band + '_ModisPredicted']))
        .addBands(date_predicted.select([band], [band + '_LSTPredicted']))
    
    var clustering = function(img) {

        // Make the training dataset.
        var training = img.select(band + '_LSTPair').sample({
            region: roi,
            scale: 30
        });

        // Instantiate the clusterer and train it.
        var clusterer = ee.Clusterer.wekaKMeans(k).train(training);

        // Cluster the input 
        var result = ee.Image(img).cluster(clusterer);
        // Range Clusters
        var group = ee.List(img.addBands(result, ['cluster']).select(band + '_LSTPair', 'cluster').reduceRegion({
            reducer: ee.Reducer.minMax().group({
                groupField: 1,
                groupName: 'cluster',
            }),
            geometry: roi,
            scale: 30
        }).get('groups'))//.aside(print,'Range Clusters '+band)


        return img.addBands(result, ['cluster'])
    }
    
    var candidates = function(name) {
        return cluster_img.select(name).neighborhoodToBands(
            ee.Kernel.square({
                radius: ee.Number(window_pixels).subtract(1).divide(2),
                units: 'pixels'
            })).toArray().rename(name + '_candidates')
    }

    var euclideanDistance = function(img, original) {

        var ed = img.expression('sqrt(pow((ltc-lt),2) + pow((abs(lgc)-abs(lg)),2))', {
            ltc: img.select(['latitude_candidates']),
            lgc: img.select(['longitude_candidates']),
            lt: original.select(['latitude']),
            lg: original.select(['longitude'])
        }).rename('distance')

        return ee.Image(1).add(ed.select('distance').divide(60)).rename('distance');
    }
    var normalizedC = function(img, original) {

        var CSpectral = (img.select('SpectralDiff_candidates').multiply(10000).add(1)).log()
        var CTemp = (img.select('TempDiff_candidates').multiply(10000).add(1)).log()
        var c = CSpectral.multiply(CTemp).multiply(original.select('distance'))

        return ee.Image(1).divide(c).rename('NormalizedC')

    }


    var cluster_img = clustering(lst_img.unmask(0));


    var allCandidatesArrays = cluster_img.addBands(candidates('cluster'))
                                          .addBands(candidates(band + '_LSTPair'))
                                          .addBands(candidates(band + '_ModisPair'))
                                          .addBands(candidates(band + '_ModisPredicted'))
                                          .addBands(candidates('latitude'))
                                          .addBands(candidates('longitude'))
                                          .addBands(candidates('TempDiff'))
                                          .addBands(candidates('SpectralDiff'))

    var mask_clusters = (allCandidatesArrays.select('cluster').eq(allCandidatesArrays.select('cluster_candidates')))
    //.and(allCandidatesArrays.select('SpectralDiff_candidates').lt(lst_img.select('SpectralDiff')))


    var candidates_img = allCandidatesArrays.select(band + '_LSTPair_candidates', 'cluster_candidates', band + '_ModisPair_candidates',
        band + '_ModisPredicted_candidates', 'latitude_candidates', 'longitude_candidates', 'TempDiff_candidates', 'SpectralDiff_candidates').arrayMask(mask_clusters)
    //Map.addLayer(candidates_img.clip(ROI), {}, '')  

    var final1 = lst_img.addBands(candidates_img).addBands(euclideanDistance(candidates_img, lst_img))
    var final = final1.addBands(normalizedC(candidates_img, final1))



    var nCandidates = candidates_img.select('cluster_candidates').arrayLengths().arrayGet([0]).rename('nCandidates')
    var ndviValues = (final.select(band + '_ModisPredicted_candidates').add(final.select(band + '_LSTPair_candidates')))
        .subtract(final.select(band + '_ModisPair_candidates')).rename(band + 'Values')

    var weight = final.select('NormalizedC').divide(final.select('NormalizedC').arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0])).rename('weight')
    var finalPredicted = (weight.select('weight').multiply(ndviValues.select(band + 'Values'))).arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0]).rename(band).addBands(nCandidates);

    var error_modis = lst_img.select(band + '_ModisPredicted').subtract(finalPredicted.select(band)).rename('error_Modis')
    var error_lst = lst_img.select(band + '_LSTPredicted').subtract(finalPredicted.select(band)).rename('error_LST')
 
    return finalPredicted.addBands(lst_img.select(band + '_ModisPredicted')).addBands(lst_img.select(band + '_LSTPredicted')).addBands(cluster_img.select('cluster')).addBands(ee.Image.pixelLonLat())
        .addBands(error_modis).addBands(error_lst)

}

//////////////////////////////////////////////////////// Analysis /////////////////////////////////////////////
var rsquared = function(original, estimated, band, roi) {

    //sse =∑ni (yi−y^i)2
    var sse = original.select(band).subtract(estimated.select(band)).pow(2)

    var Sum_sse = sse.reduceRegion(ee.Reducer.sum(), roi, 30, original.select(band).projection()).get(band)
    var mean = original.reduceRegion(ee.Reducer.mean(), roi, 30, original.select(band).projection()).get(band)
    var ssto = original.select(band).subtract(ee.Number(mean)).pow(2);
    var Sum_ssto = ssto.reduceRegion(ee.Reducer.sum(), roi, 30, original.select(band).projection()).get(band)

    // Rsquare = 1 - (first sum/second sum)
    return ee.Number(1).subtract(ee.Number(Sum_sse).divide(ee.Number(Sum_ssto)))

}

var histograms = function(date_pair1, date_predicted, roi, roi_name, flag) {

    if (flag) {

        var w = [3, 5, 7, 9];
        var k = [4, 6, 8];

        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {

                var options = {
                    title: 'NDVI histogram in ' + roi_name + ' W' + nameW + 'K' + nameK,
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
                var pre = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted)

                var histogram = ui.Chart.image.histogram(pre.updateMask(pre.neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).select(['NDVI', 'NDVI_LSTPredicted']), roi, 30)
                    .setSeriesNames(['NDVI_Predicted', 'NDVI_LST'])
                    .setOptions(options);

                // Display the histogram.
                print(histogram);
            })

            var k4 = predict(w[indexW], 4, 'NDVI', roi, date_pair1, date_predicted).select(['NDVI'], ['NDVI_K4']).updateMask(date_predicted.select('NDVI'))
            var k6 = predict(w[indexW], 6, 'NDVI', roi, date_pair1, date_predicted).select(['NDVI'], ['NDVI_K6']).updateMask(date_predicted.select('NDVI'))
            var k8 = predict(w[indexW], 8, 'NDVI', roi, date_pair1, date_predicted).select(['NDVI'], ['NDVI_K8']).addBands(k4).addBands(k6).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).updateMask(k4.neq(0))

            var histogramW = ui.Chart.image.histogram(k8.select(['NDVI_K4', 'NDVI_K6', 'NDVI_K8']), roi, 30)
                .setSeriesNames(['NDVI_K4', 'NDVI_K6', 'NDVI_K8'])
                .setOptions({
                    title: 'NDVI histogram in ' + roi_name + ' W' + nameW,
                    fontSize: 10,
                    hAxis: {
                        title: 'NDVI'
                    },
                    vAxis: {
                        title: 'count of NDVI'
                    }
                });

            // Display the histogram.
            print(histogramW);

        })

        var image = predict(3, 4, 'NDVI', roi, date_pair1, date_predicted).updateMask(date_predicted.select('NDVI'))

        var histogramOriginals = ui.Chart.image.histogram(image.updateMask(image.neq(0)).select(['NDVI_LSTPredicted', 'NDVI_ModisPredicted']), roi, 30)
            .setSeriesNames(['NDVI_LSTPredicted', 'NDVI_ModisPredicted'])
            .setOptions({
                title: 'Original NDVI Histogram in ' + roi_name,
                fontSize: 10,
                hAxis: {
                    title: 'NDVI'
                },
                vAxis: {
                    title: 'count of NDVI'
                },
                series: {
                    0: {
                        color: '#1c9099'
                    },
                    1: {
                        color: '#756bb1'
                    }
                }
            });
        print(histogramOriginals)

    } else {

        var w = [3, 5, 7, 9];
        var k = [4, 6, 8];

        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {

                var options = {
                    title: 'NDVI histogram in ' + roi_name + ' W' + nameW + 'K' + nameK,
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

                var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_predicted);
                var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_predicted);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                var error_LST = date_predicted.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')
                var pre = ndvi.addBands(date_predicted.select(['NDVI'], ['NDVI_LSTPredicted']))

                var histogram = ui.Chart.image.histogram(pre.updateMask(pre.neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).select(['NDVI', 'NDVI_LSTPredicted']), roi, 30)
                    .setSeriesNames(['NDVI_Predicted', 'NDVI_LST'])
                    .setOptions(options);

                // Display the histogram.
                print(histogram);
            })

            // K=4
            var imageNIR4 = predict(w[indexW], 4, 'nir', roi, date_pair1, date_predicted);
            var imageRED4 = predict(w[indexW], 4, 'red', roi, date_pair1, date_predicted);
            var k4 = imageNIR4.select('nir').addBands(imageRED4.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI_K4').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
            // K=6
            var imageNIR6 = predict(w[indexW], 6, 'nir', roi, date_pair1, date_predicted);
            var imageRED6 = predict(w[indexW], 6, 'red', roi, date_pair1, date_predicted);
            var k6 = imageNIR6.select('nir').addBands(imageRED6.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI_K6').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))

            // K=8
            var imageNIR8 = predict(w[indexW], 8, 'nir', roi, date_pair1, date_predicted);
            var imageRED8 = predict(w[indexW], 8, 'red', roi, date_pair1, date_predicted);
            var k8 = imageNIR8.select('nir').addBands(imageRED8.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI_K8').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))


            var predicted = k8.addBands(k4).addBands(k6).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).updateMask(k4.neq(0))


            var histogramW = ui.Chart.image.histogram(predicted.select(['NDVI_K4', 'NDVI_K6', 'NDVI_K8']), roi, 30)
                .setSeriesNames(['NDVI_K4', 'NDVI_K6', 'NDVI_K8'])
                .setOptions({
                    title: 'NDVI histogram in ' + roi_name + ' W' + nameW,
                    fontSize: 10,
                    hAxis: {
                        title: 'NDVI'
                    },
                    vAxis: {
                        title: 'count of NDVI'
                    }
                });

            // Display the histogram.
            print(histogramW);

        })

        var image = predict(3, 4, 'NDVI', roi, date_pair1, date_predicted).updateMask(date_predicted.select('NDVI'))

        var histogramOriginals = ui.Chart.image.histogram(image.updateMask(image.neq(0)).select(['NDVI_LSTPredicted', 'NDVI_ModisPredicted']), roi, 30)
            .setSeriesNames(['NDVI_LSTPredicted', 'NDVI_ModisPredicted'])
            .setOptions({
                title: 'Original NDVI Histogram in ' + roi_name,
                fontSize: 10,
                hAxis: {
                    title: 'NDVI'
                },
                vAxis: {
                    title: 'count of NDVI'
                },
                series: {
                    0: {
                        color: '#1c9099'
                    },
                    1: {
                        color: '#756bb1'
                    }
                }
            });
        print(histogramOriginals)

    }

}

var scatters = function(date_pair1, date_predicted, roi, roi_name, flag) {

    var w = [3, 5, 7, 9];
    var k = [4, 6, 8];

    if (flag) {
        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {

                var predicted = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted)
                predicted = predicted.updateMask(predicted.select('NDVI').neq(0))

                var original = date_predicted.updateMask(predicted.select('NDVI').unmask(-2).neq(-2))
                    .addBands(predicted.select(['NDVI'], ['NDVI_Predicted']).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)))
                    .reduceRegion(ee.Reducer.toList(), roi, 30)

                var x = original.get('NDVI');
                var y = original.get('NDVI_Predicted');
                var chart = ui.Chart.array.values(y, 0, x);

                chart = chart.setSeriesNames(['NDVI'])
                chart = chart.setOptions({
                    title: 'NDVI ScatterPlot in ' + roi_name + ' W' + nameW + 'K' + nameK,
                    hAxis: {
                        title: "LST Original NDVI",
                        viewWindow: {
                            min: 0,
                            max: 1
                        }
                    },
                    vAxis: {
                        title: 'Predicted NDVI W' + nameW + 'K' + nameK,
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

            })
        })


    } else {

        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {

                var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_predicted);
                var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_predicted);
                var predicted = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))

                predicted = predicted.updateMask(predicted.select('NDVI').neq(0))

                var original = date_predicted.updateMask(predicted.select('NDVI').unmask(-2).neq(-2))
                    .addBands(predicted.select(['NDVI'], ['NDVI_Predicted']).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)))
                    .reduceRegion(ee.Reducer.toList(), roi, 30)

                var x = original.get('NDVI');
                var y = original.get('NDVI_Predicted');
                var chart = ui.Chart.array.values(y, 0, x);

                chart = chart.setSeriesNames(['NDVI'])
                chart = chart.setOptions({
                    title: 'NDVI ScatterPlot in ' + roi_name + ' W' + nameW + 'K' + nameK,
                    hAxis: {
                        title: "LST Original NDVI",
                        viewWindow: {
                            min: 0,
                            max: 1
                        }
                    },
                    vAxis: {
                        title: 'Predicted NDVI W' + nameW + 'K' + nameK,
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

            })
        })
    }

}

var clusters = function(img, band, k, roi) {

    // Make the training dataset.
    var training = img.select(band).sample({
        region: roi,
        scale: 30
    });

    // Instantiate the clusterer and train it.
    var clusterer = ee.Clusterer.wekaKMeans(k).train(training);


    // Cluster the input 
    var result = ee.Image(img).cluster(clusterer);

    // Range Clusters
    var group = ee.List(img.addBands(result, ['cluster']).select(band, 'cluster').reduceRegion({
        reducer: ee.Reducer.minMax().group({
            groupField: 1,
            groupName: 'cluster',
        }),
        geometry: roi,
        scale: 30
    }).get('groups')).aside(print, 'Range Clusters ' + band)

    return img.addBands(result, ['cluster'])
}

///////////////////////////// APP /////////////////////////////////

var app = {};

/** Creates the UI panels. */
app.createPanels = function() {

    app.intro = {
        panel: ui.Panel([
            ui.Label({
                value: 'STARFM 1Pair',
                style: {
                    fontWeight: 'bold',
                    fontSize: '24px',
                    margin: '10px 5px',
                    shown: true
                }
            }),
            ui.Label('This app allows you to apply STARFM using 1Pair Landsat-Modis')
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

    /* The panel for the visualization section with corresponding widgets. */
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

    /* The collection filter controls. */
    app.filters = {
        startDate: ui.Textbox('YYYY-MM-DD', '2017-01-01'),
        endDate: ui.Textbox('YYYY-MM-DD', '2017-12-31'),
        applyButton: ui.Button('Apply filters', app.applyDates),
        loadingLabel: ui.Label({
            value: 'Loading...',
            style: {
                stretch: 'vertical',
                color: 'gray',
                shown: false
            }
        })
    };

    /* The panel for the filter control widgets. */
    app.filters.panel = ui.Panel({
        widgets: [
            ui.Label('1) Select Dates of Interest', {
                fontWeight: 'bold'
            }),
            ui.Label('Start date', app.HELPER_TEXT_STYLE), app.filters.startDate,
            ui.Label('End date', app.HELPER_TEXT_STYLE), app.filters.endDate,
            ui.Panel([
                app.filters.applyButton,
                app.filters.loadingLabel
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

    /* The image picker section. */
    app.picker = {
        // Create a select with a function that reacts to the "change" event.
        predictedDate: {
            label: ui.Label(),

            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.predictedDate.label.setValue('Cloud Cover in ROI: ' + result.toFixed(2))
                    })
                }
            })
        },
        pairDate: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.pairDate.label.setValue('Cloud Cover in ROI: ' + result.toFixed(2))
                    })
                }
            })
        }
    };

    /* The panel for the picker section with corresponding widgets. */
    app.picker.panel = ui.Panel({
        widgets: [
            ui.Label('2) Select Images', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Pair Date', app.HELPER_TEXT_STYLE), app.picker.pairDate.select,
                app.picker.pairDate.label
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Predicted Date', app.HELPER_TEXT_STYLE), app.picker.predictedDate.select,
                app.picker.predictedDate.label
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

    app.KW = {
        applyNDVI: ui.Checkbox({
            label: 'Apply NDVI at first',
            value: true
        }),
        K: {
            label: ui.Label(),
            // Create a select with a function that reacts to the "change" event.
            select: ui.Select({
                items: Object.keys(app.K_OPTIONS),
                onChange: function() {
                    // Update the label's value with the select's description.
                    var option = app.K_OPTIONS[app.KW.K.select.getValue()];
                    app.KW.K.label.setValue(option.description);
                }
            })
        },
        W: {
            label: ui.Label(),
            // Create a select with a function that reacts to the "change" event.
            select: ui.Select({
                items: Object.keys(app.W_OPTIONS),
                onChange: function() {
                    // Update the label's value with the select's description.
                    var option = app.W_OPTIONS[app.KW.W.select.getValue()];
                    app.KW.W.label.setValue(option.description);

                }
            })
        }
    };

    /* The panel for the visualization section with corresponding widgets. */
    app.KW.panel = ui.Panel({
        widgets: [
            ui.Label('3) Select K; W, NDVI', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('K', app.HELPER_TEXT_STYLE), app.KW.K.select,
                ui.Label('W', app.HELPER_TEXT_STYLE), app.KW.W.select,
            ], ui.Panel.Layout.flow('horizontal')),
            app.KW.applyNDVI
        ],
        style: app.SECTION_STYLE
    })


    // Default the select to the first value.
    app.KW.K.select.setValue(app.KW.K.select.items().get(0));
    app.KW.W.select.setValue(app.KW.W.select.items().get(0));

    app.show = {
        buttonPairs: ui.Button({
            label: 'Show Pair LST-Modis',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapPairs()
            }
        }),
        buttonPairPredicted: ui.Button({
            label: 'Show Pair LST-Predicted',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapPairPredicted()
            }
        }),
        buttonOriginalLST: ui.Button({
            label: 'Show Original LST',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalLST()
            }
        }),
        buttonPredicted: ui.Button({
            label: 'Show Predicted',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapPredicted()
            }
        }),
        buttonClusters: ui.Button({
            label: 'Show Clusters',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapClusters()
            }
        }),
        buttonErrors: ui.Button({
            label: 'Show Errors',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapErrors()
            }
        }),
        buttonOriginalPredicted: ui.Button({
            label: 'Show Predicted - Original LST',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapCompare()
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
        buttonOriginalModis: ui.Button({
            label: 'Show Original Modis',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalModis()
            }
        }),
        buttonInspector: ui.Button({
            label: 'Inspector NDVI',
            // React to the button's click event.
            onClick: function() {
                app.inspectorNDVI()
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

    /* The panel for the export section with corresponding widgets. */
    app.show.panel = ui.Panel({
        widgets: [
            ui.Label('4) Visualization', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Compare', app.HELPER_TEXT_STYLE),
                app.show.buttonPairs,
                app.show.buttonOriginalPredicted,
                app.show.buttonPairPredicted
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Individual', app.HELPER_TEXT_STYLE),
                app.show.buttonPredicted,
                app.show.buttonOriginalLST,
                app.show.buttonOriginalModis
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Other', app.HELPER_TEXT_STYLE),
                app.show.buttonClusters,
                app.show.buttonErrors
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Analysis', app.HELPER_TEXT_STYLE),
                app.show.buttonStatsRegion,
                app.show.buttonMetrics,
                app.show.buttonHistograms,
                app.show.buttonScatters,
                app.show.buttonInspector
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

    /* The export section. */
    app.export = {
        buttonPairs: ui.Button({
            label: 'Export Pair Images',
            // React to the button's click event.
            onClick: function() {
                app.exportPairs();
            }
        }),
        buttonOriginals: ui.Button({
            label: 'Export Original Images',
            // React to the button's click event.
            onClick: function() {
                app.exportOriginals();
            }
        }),
        buttonPredicted: ui.Button({
            label: 'Export Predicted Image',
            // React to the button's click event.
            onClick: function() {
                app.exportPredictedImage();
            }
        })
    };

    /* The panel for the export section with corresponding widgets. */
    app.export.panel = ui.Panel({
        widgets: [
            ui.Label('5) Start an export', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                app.export.buttonPredicted,
                app.export.buttonOriginals,
                app.export.buttonPairs
            ], ui.Panel.Layout.flow('horizontal'))
        ],
        style: app.SECTION_STYLE
    });

}

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
            app.KW.K.select,
            app.KW.W.select,
            app.filters.applyButton,
            app.picker.predictedDate.select,
            app.picker.pairDate.select
        ];
        loadDependentWidgets.forEach(function(widget) {
            widget.setDisabled(enabled);
        });
    };

    /** Applies the selection filters currently selected in the UI. */
    app.applyDates = function() {
        app.setLoadingMode(true);
        var ROI = app.ROI_OPTIONS[app.roi.select.getValue()].roi

        var filtered = ee.ImageCollection(L8)
            .filterBounds(ROI)

        // Set filter variables.
        var start = app.filters.startDate.getValue();
        if (start) start = ee.Date(start);
        var end = app.filters.endDate.getValue();
        if (end) end = ee.Date(end);
        if (start) filtered = filtered.filterDate(start, end)
            .map(maskCloudsLandsat)
            .select(['B5', 'B4'], ['nir', 'red'])
            .map(addNDVI)
            .map(addTime)
            .map(function(img) {
                var clipped = img.clip(ROI);
                var img_nulls = clipped.unmask(-2);
                var total_pixels = ee.Dictionary(['TOTAL_PIXELS', clipped.unmask(-2).reduceRegion(ee.Reducer.count(), ROI, 30).get('NDVI')]);
                var null_pixels = ee.Dictionary(['NULL_PIXELS', img_nulls.updateMask(img_nulls.eq(-2)).reduceRegion(ee.Reducer.count(), ROI, 30).get('NDVI')]);
                var quality_pixels = ee.Dictionary(['QUALITY_PIXELS', clipped.reduceRegion(ee.Reducer.count(), ROI, 30).get('NDVI')]);
                var cloud_cover_region = ee.Dictionary(['CLOUD_COVER_REGION', ee.Number(null_pixels.get('NULL_PIXELS')).divide(total_pixels.get('TOTAL_PIXELS')).multiply(100)]);
                var dict = total_pixels.combine(null_pixels).combine(quality_pixels).combine(cloud_cover_region);
                var image = clipped.set(dict);

                return image.addBands(ee.Image.constant(null_pixels.get('NULL_PIXELS')).rename('Null_Pixels'))
                    .addBands(ee.Image.constant(cloud_cover_region.get('CLOUD_COVER_REGION')).rename('Cloud_Cover_Region'))
                    .addBands(ee.Image.constant(quality_pixels.get('QUALITY_PIXELS')).rename('Quality_Pixels'))
                    .addBands(ee.Image.pixelLonLat())
            })

        // Get the list of computed ids.
        var computedIds = filtered
            .limit(22)
            .reduceColumns(ee.Reducer.toList(), ['system:index'])
            .get('list');

        computedIds.evaluate(function(ids) {
            // Update the image picker with the given list of ids.
            app.setLoadingMode(false);
            app.picker.predictedDate.select.items().reset(ids);
            app.picker.pairDate.select.items().reset(ids);
            // Default the image picker to the first id.
            app.picker.predictedDate.select.setValue(app.picker.predictedDate.select.items().get(0));
            app.picker.pairDate.select.setValue(app.picker.pairDate.select.items().get(0));
        });

        collection = filtered;
        print(collection)

    };

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

    app.inspectorNDVI = function() {

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

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Original LST Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')
        var predictedDate = ee.Image(app.COLLECTION_ID + '/' + app.picker.predictedDate.select.getValue());
        if (predictedDate) {
            // If an image id is found, create an image.
            var image = ee.ImageCollection(predictedDate).map(maskClouds).select(['B5', 'B4'], ['nir', 'red']).map(addNDVI)
            map.addLayer(ee.Image(image.first()).select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            map.addLayer(ee.Image(image.first()).select('NDVI').unmask(-1).eq(-1).updateMask(ee.Image(image.first()).select('NDVI').unmask(-1).eq(-1)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(ee.Image(image.first()).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original')

        }

        var inspector = ui.Panel({
            layout: ui.Panel.Layout.flow('horizontal'),
            style: {
                position: 'bottom-right'
            }
        });

        var intro = ui.Panel([
            ui.Label('Click a point on the map to inspect.')
        ]);

        inspector.add(intro);

        var ndvi = ee.ImageCollection(app.COLLECTION_ID)
            .filterBounds(app.ROI_OPTIONS[app.roi.select.getValue()].roi)
            .filterDate(app.filters.startDate.getValue(), app.filters.endDate.getValue())
            .map(maskClouds).select(['B5', 'B4'], ['nir', 'red']).map(addNDVI)

        // Create panels to hold lon/lat values.
        var lon = ui.Label();
        var lat = ui.Label();
        inspector.add(ui.Panel([lon, lat], ui.Panel.Layout.flow('vertical')));
        // Register a callback on the default map to be invoked when the map is clicked.
        map.onClick(function(coords) {
            intro.clear();
            // Update the lon/lat panel with values from the click event.
            lon.setValue('lon: ' + coords.lon.toFixed(2)),
                lat.setValue('lat: ' + coords.lat.toFixed(2));

            // Add a red dot for the point clicked on.
            var point = ee.Geometry.Point(coords.lon, coords.lat);
            var dot = ui.Map.Layer(point, {
                color: 'FF0000'
            });
            //map.layers().set(1, dot);
            // Create an NDVI chart.
            var ndviChart = ui.Chart.image.series(ndvi.select('NDVI'), point, ee.Reducer.mean(), 30);
            ndviChart.setOptions({
                title: 'NDVI Over Time',
                vAxis: {
                    title: 'NDVI'
                },
                hAxis: {
                    title: 'date',
                    format: 'MM-yy',
                    gridlines: {
                        count: 7
                    }
                },
            });
            inspector.widgets().set(2, ndviChart);

        });
        map.style().set('cursor', 'crosshair');
        map.add(inspector);
    }

    app.refreshMapCompare = function() {

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
        var originalDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first());

        var NAMES = [
            'Original L8 Image',
            'PredictedImage'
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

        if (originalDate) {
            // If an image id is found, create an image.
            maps[0].addLayer(originalDate.select('NDVI'), ndviParams, 'OriginalImage');
            maps[0].addLayer(originalDate.select('NDVI').unmask(-2).eq(-2).updateMask(originalDate.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');

        }

        if (pairDate && originalDate) {

            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var predicted = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', roi,
                    pairDate, originalDate);

                maps[1].addLayer(predicted.select('NDVI'), ndviParams, 'PredictedImage');
                maps[1].addLayer(predicted.select('NDVI').eq(0).updateMask(predicted.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)
                // Add a label to the panel.
                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = predicted.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var originalPoint = originalDate.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');
                    var computedValueOriginal = originalPoint.get('NDVI')

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        computedValueOriginal.evaluate(function(result) {
                            // Add a label with the results from the server.
                            inspector.add(ui.Label({
                                value: 'NDVI original: ' + result.toFixed(2),
                                style: {
                                    stretch: 'horizontal'
                                }
                            }));
                        })
                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI Predicted: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));
                    });
                });
            } else {

                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', roi,
                    pairDate, originalDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', roi,
                    pairDate, originalDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI')
                    .addBands(originalDate.select(['NDVI'], ['NDVI_LSTPredicted']));
                maps[1].addLayer(ndvi.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                maps[1].addLayer(ndvi.select('NDVI').eq(0).updateMask(ndvi.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)
                // Add a label to the panel.
                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = ndvi.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var originalPoint = originalDate.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');
                    var computedValueOriginal = originalPoint.get('NDVI')

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        computedValueOriginal.evaluate(function(result) {
                            // Add a label with the results from the server.
                            inspector.add(ui.Label({
                                value: 'NDVI original: ' + result.toFixed(2),
                                style: {
                                    stretch: 'horizontal'
                                }
                            }));
                        })
                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI Predicted: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });
            }
        }
        ui.root.widgets().set(0, mapGrid);

    };

    app.refreshMapPairPredicted = function() {

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
        var originalDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first());
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first());

        var NAMES = [
            'Pair L8 Image',
            'PredictedImage'
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

        if (pairDate) {
            // If an image id is found, create an image.
            maps[0].addLayer(pairDate.select('NDVI'), ndviParams, 'pairImage');
            maps[0].addLayer(pairDate.select('NDVI').unmask(-2).eq(-2).updateMask(pairDate.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');

        }

        if (pairDate && originalDate) {

            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var predicted = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', roi,
                    pairDate, originalDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                maps[1].addLayer(predicted.select('NDVI'), ndviParams, 'PredictedImage');
                maps[1].addLayer(predicted.select('NDVI').eq(0).updateMask(predicted.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)


                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = predicted.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var originalPoint = predictedDate.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');
                    var computedValueOriginal = originalPoint.get('NDVI')

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        computedValueOriginal.evaluate(function(result) {
                            // Add a label with the results from the server.
                            inspector.add(ui.Label({
                                value: 'NDVI pair: ' + result.toFixed(2),
                                style: {
                                    stretch: 'horizontal'
                                }
                            }));
                        })
                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI predicted: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });


            } else {

                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', roi,
                    pairDate, originalDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', roi,
                    pairDate, originalDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI')
                    .addBands(originalDate.select(['NDVI'], ['NDVI_LSTPredicted']));
                maps[1].addLayer(ndvi.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                maps[1].addLayer(ndvi.select('NDVI').eq(0).updateMask(ndvi.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)


                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = ndvi.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var originalPoint = predictedDate.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');
                    var computedValueOriginal = originalPoint.get('NDVI')

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        computedValueOriginal.evaluate(function(result) {
                            // Add a label with the results from the server.
                            inspector.add(ui.Label({
                                value: 'NDVI pair: ' + result.toFixed(2),
                                style: {
                                    stretch: 'horizontal'
                                }
                            }));
                        })
                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI predicted: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });

            }

        }
        ui.root.widgets().set(0, mapGrid);

    }

    app.refreshMapOriginalLST = function() {

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

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Original LST Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first());
        if (predictedDate) {
            // If an image id is found, create an image.
            map.addLayer(predictedDate.select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            map.addLayer(predictedDate.select('NDVI').unmask(-2).eq(-2).updateMask(predictedDate.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(predictedDate.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original')


            var inspector = ui.Panel({
                layout: ui.Panel.Layout.flow('horizontal'),
                style: {
                    position: 'bottom-right'
                }
            });
            inspector.add(ui.Label('Click to get NDVI'));

            // Add the panel to the default map.
            map.add(inspector);

            // Set the default map's cursor to a "crosshair".
            map.style().set('cursor', 'crosshair');

            // Register a callback on the map to be invoked when the map is clicked.
            map.onClick(function(coords) {
                // Clear the panel and show a loading message.
                inspector.clear();
                inspector.style().set('shown', true);
                inspector.add(ui.Label('Loading...', {
                    color: 'gray'
                }));

                // Compute the mean NDVI; a potentially long-running server operation.
                var point = ee.Geometry.Point(coords.lon, coords.lat);
                var sampledPoint = predictedDate.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                var computedValue = sampledPoint.get('NDVI');

                // Request the value from the server and use the results in a function.
                computedValue.evaluate(function(result) {
                    inspector.clear();

                    // Add a label with the results from the server.
                    inspector.add(ui.Label({
                        value: 'NDVI: ' + result.toFixed(2),
                        style: {
                            stretch: 'vertical'
                        }
                    }));

                });
            });

        }
    };

    app.refreshMapOriginalModis = function() {

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
        map.add(ui.Label('Predicted Modis Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')

        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).map(addTime).first())

        var modisPredicted = ee.Image(modis
            .filterBounds(roi)
            .map(addTime)
            .filterMetadata('Date', 'equals', predictedDate.get('Date'))
            .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
            .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
            crs: predictedDate.select('nir').projection()
        })

        map.addLayer(modisPredicted.select('NDVI').clip(roi), ndviParams, 'OriginalModisPair');
        print(modisPredicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original Modis')

        var inspector = ui.Panel({
            layout: ui.Panel.Layout.flow('horizontal'),
            style: {
                position: 'bottom-right'
            }
        });
        inspector.add(ui.Label('Click to get NDVI'));

        // Add the panel to the default map.
        map.add(inspector);

        // Set the default map's cursor to a "crosshair".
        map.style().set('cursor', 'crosshair');

        // Register a callback on the map to be invoked when the map is clicked.
        map.onClick(function(coords) {
            // Clear the panel and show a loading message.
            inspector.clear();
            inspector.style().set('shown', true);
            inspector.add(ui.Label('Loading...', {
                color: 'gray'
            }));

            // Compute the mean NDVI; a potentially long-running server operation.
            var point = ee.Geometry.Point(coords.lon, coords.lat);
            var sampledPoint = modisPredicted.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
            var computedValue = sampledPoint.get('NDVI');

            // Request the value from the server and use the results in a function.
            computedValue.evaluate(function(result) {
                inspector.clear();

                // Add a label with the results from the server.
                inspector.add(ui.Label({
                    value: 'NDVI: ' + result.toFixed(2),
                    style: {
                        stretch: 'vertical'
                    }
                }));

            });
        });
    }

    app.refreshMapPredicted = function() {

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

        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Predicted Image'))
        map.add(legend)
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI')
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())


        if (pairDate && predictedDate) {
            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var image = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                map.addLayer(image.select('NDVI').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), ndviParams, 'PredictedImage');
                map.addLayer(image.select('NDVI').eq(0).updateMask(image.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                print(image.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                map.add(inspector);

                // Set the default map's cursor to a "crosshair".
                map.style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                map.onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = image.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });

            } else {
                // If an image id is found, create an image.
                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                map.addLayer(ndvi.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                map.addLayer(ndvi.select('NDVI').eq(0).updateMask(ndvi.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                print(ndvi.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                map.add(inspector);

                // Set the default map's cursor to a "crosshair".
                map.style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                map.onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = ndvi.select('NDVI').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('NDVI');

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'NDVI: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });

            }
        }

    };

    app.refreshMapPairs = function() {

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
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).map(addTime).first())

        var NAMES = [
            'Pair LST',
            'Pair Modis'
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

        if (pairDate) {
            maps[0].addLayer(pairDate.select('NDVI').clip(roi), ndviParams, 'PairLST');
            maps[0].addLayer(pairDate.select('NDVI').unmask(-2).eq(-2).updateMask(pairDate.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(pairDate.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair LST')

            var modisPair = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
                .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
                crs: pairDate.select('NDVI').projection()
            })

            print(modisPair.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair Modis')
            maps[1].addLayer(modisPair.select('NDVI').clip(roi), ndviParams, 'PairModis');
            maps[1].add(legend);
        }

        ui.root.widgets().set(0, mapGrid);


    };

    app.refreshMapClusters = function() {
        Map.clear()
        var map = ui.Map()
        ui.root.widgets().set(0, map);
        map.add(ui.Label('Clusters - Pair Image'))
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var k = app.K_OPTIONS[app.KW.K.select.getValue()].k
        map.centerObject(app.ROI_OPTIONS[app.roi.select.getValue()].roi, 12);
        map.addLayer(app.ROI_OPTIONS[app.roi.select.getValue()].roi, {}, 'ROI')

        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())

        if (app.KW.applyNDVI.getValue()) {
            if (pairDate && predictedDate) {

                map.addLayer(clusters(pairDate, 'NDVI', k, roi).select('cluster').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi).randomVisualizer(), {}, 'Clusters-PredictedImage');
            }
        } else {
            print('Impossible to Show Clusters. Select apply NDVI from the beginning')
        }
    };

    app.refreshMapErrors = function() {

        // Create the panel for the legend items.
        var legend = ui.Panel({
            style: {
                // width: '200px',
                position: 'bottom-left',
                padding: '8px 15px'
            }
        });

        legend.add(addColor('#fc8d59', '-0.2')).add(addColor('#ffffbf', '')).add(addColor('#91cf60', '0.2'))

        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var originalDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())

        var NAMES = [
            'Error Modis',
            'Error LST'
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

        if (pairDate && originalDate) {

            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var predicted = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, originalDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                maps[0].addLayer(predicted.select('error_Modis').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), errorParams, 'Error_Modis');
                maps[1].addLayer(predicted.select('error_LST').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), errorParams, 'Error_LST');
                maps[1].add(legend)
                print(predicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = predicted.select('error_LST').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('error_LST');

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'Error_NDVI_LST: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });

            } else {

                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', roi,
                    pairDate, originalDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', roi,
                    pairDate, originalDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                var error_LST = originalDate.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')

                var modisPredicted = ee.Image(modis
                    .filterBounds(roi)
                    .map(addTime)
                    .filterMetadata('Date', 'equals', originalDate.get('Date'))
                    .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
                    .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
                    crs: originalDate.select('nir').projection()
                })

                var error_Modis = modisPredicted.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_Modis')
                maps[1].addLayer(error_LST.select('error_LST').clip(roi), errorParams, 'Error_LST');
                maps[0].addLayer(error_Modis.select('error_Modis').clip(roi), errorParams, 'Error_Modis');
                maps[1].add(legend)
                print(error_LST.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean error_LST')
                print(error_Modis.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean error_Modis')

                var inspector = ui.Panel({
                    layout: ui.Panel.Layout.flow('horizontal'),
                    style: {
                        position: 'bottom-right'
                    }
                });
                inspector.add(ui.Label('Click to get NDVI'));

                // Add the panel to the default map.
                maps[1].add(inspector);

                // Set the default map's cursor to a "crosshair".
                maps[1].style().set('cursor', 'crosshair');

                // Register a callback on the map to be invoked when the map is clicked.
                maps[1].onClick(function(coords) {
                    // Clear the panel and show a loading message.
                    inspector.clear();
                    inspector.style().set('shown', true);
                    inspector.add(ui.Label('Loading...', {
                        color: 'gray'
                    }));

                    // Compute the mean NDVI; a potentially long-running server operation.
                    var point = ee.Geometry.Point(coords.lon, coords.lat);
                    var sampledPoint = error_LST.select('error_LST').reduceRegion(ee.Reducer.first(), point, 30);
                    var computedValue = sampledPoint.get('error_LST');

                    // Request the value from the server and use the results in a function.
                    computedValue.evaluate(function(result) {
                        inspector.clear();

                        // Add a label with the results from the server.
                        inspector.add(ui.Label({
                            value: 'Error_NDVI: ' + result.toFixed(2),
                            style: {
                                stretch: 'vertical'
                            }
                        }));

                    });
                });

            }
            ui.root.widgets().set(0, mapGrid);
        };
    }

    app.showHistograms = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())


        print('--------- Histograms--------------')
        histograms(date_pair1, date_predicted, roi, roi_name, app.KW.applyNDVI.getValue());
    }

    app.showScatters = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())

        print('--------- Scatter-----------------')
        scatters(date_pair1, date_predicted, roi, roi_name, app.KW.applyNDVI.getValue());

    }

    app.metrics = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())

        app.exportMetrics();
        var w = [3, 5, 7, 9];
        var k = [4, 6, 8];
        var dici = ee.Dictionary();

        if (app.KW.applyNDVI.getValue()) {

            print('----------------R Squared NDVI-----------------')
            w.forEach(function(nameW, indexW) {
                k.forEach(function(nameK, indexK) {
                    print('---------K=' + nameK + '------------')
                    var image = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted)
                    image = image.updateMask(image.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))

                    print(rsquared(date_predicted.updateMask(image.select('NDVI').neq(0)), image, 'NDVI', roi), 'R squared W' + nameW + 'K' + nameK)

                    // Original - Predicted
                    var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
                    var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
                    var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
                    var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
                    var corr = image.select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
                    var dict = dici.combine(mae).combine(rmse).combine(error_dict).combine(corr)
                    print(dict, 'Stats W' + nameW + 'K' + nameK)
                    print(predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted).updateMask(image.select('NDVI').neq(0)).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

                    print(corr, 'pearsons Correlation');

                })
            })
        } else {
            print('----------- R Squared Reflectance--------------')
            w.forEach(function(nameW, indexW) {
                k.forEach(function(nameK, indexK) {

                    print('---------K=' + nameK + '------------')
                    var imageNIR = predict(w[indexW], k[indexK], 'nir', roi,
                        date_pair1, date_predicted);
                    var imageRED = predict(w[indexW], k[indexK], 'red', roi,
                        date_pair1, date_predicted);
                    var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                    var error_LST = date_predicted.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')

                    var image = ndvi.addBands(error_LST).updateMask(ndvi.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
                        .addBands(date_predicted.select(['NDVI'], ['NDVI_LSTPredicted']))

                    print(rsquared(date_predicted.updateMask(image.select('NDVI').neq(0)), image, 'NDVI', roi), 'R squared W' + nameW + 'K' + nameK)

                    // Original - Predicted
                    var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
                    var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
                    var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
                    var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])

                    var dict = mae.combine(rmse).combine(error_dict);
                    print(dict, 'Stats W' + nameW + 'K' + nameK)
                    print(predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted).updateMask(image.select('NDVI').neq(0)).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

                    print(image.addBands(date_predicted.select(['NDVI'], ['NDVI_LSTPredicted'])).select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30), 'pearsons Correlation');

                })
            })

        }

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

    app.exportMetrics = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_abre = app.ROI_OPTIONS[app.roi.select.getValue()].abre
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
        
        var w = [3, 5, 7, 9];
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
        
        var calc_Metrics = function(date_predicted, image, nameW, nameK) {
            var r = ee.Number(rsquared(date_predicted.updateMask(image.select('NDVI').neq(0)), image, 'NDVI', roi))
  
            // Original - Predicted
            var r_dict = ee.Dictionary(['R^2', r])
            var w_dict = ee.Dictionary(['W', nameW])
            var k_dict = ee.Dictionary(['K', nameK])
            var signal_dict = ee.Dictionary(['signal', 1])
            var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
            var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
            var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
            var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])
            var corr = image.select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30)
            var dict = w_dict.combine(k_dict).combine(r_dict).combine(mae).combine(rmse).combine(error_dict).combine(corr)
          
            return ee.Feature(null, dict);
        }
        
        var makePrediction = function(indexW, indexK, nameW, nameK, date_pair1, date_predicted) {
          var feats = ee.List([]);
          // NDVI
          var imageNDVI = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_predicted)
          
          var image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
          
          var ft = calc_Metrics(date_predicted, image, nameW, nameK)
          feats = feats.cat(ee.List([ft]))
          
          //return ft;

          // Reflectance
          var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_predicted);
          var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_predicted);
          var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
          var error_LST = date_predicted.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')
          
          image = ndvi.addBands(error_LST).updateMask(ndvi.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
              .addBands(date_predicted.select(['NDVI'], ['NDVI_LSTPredicted']))
          
          ft = calc_Metrics(date_predicted, image, nameW, nameK)
          feats = feats.cat(ee.List([ft]))
          
          return feats;
        }
        
        var createBuffer = function(feature){ return feature.buffer(20); }
        
        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {
              
              var watermask = date_pair1.clip(osmwater).unmask().not();
              var roadmask = date_pair1.clip(osmroads.map(createBuffer)).unmask().not();
              var buildingmask = date_pair1.clip(osmbuildings).unmask().not();
              
              // NO MASK
              var ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1, date_predicted)
              features = features.cat(ft)
              
              // WATER MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask), date_predicted)
              featureswatermasked = featureswatermasked.cat(ft)
              
              // ROAD MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(roadmask), date_predicted)
              featuresroadsmasked = featuresroadsmasked.cat(ft)
              
              // BUILDING MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(buildingmask), date_predicted)
              featuresbuildingsmasked = featuresbuildingsmasked.cat(ft)
              
              // WR MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(roadmask), date_predicted)
              featureswrmasked = featureswrmasked.cat(ft)
              
              // WB MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(buildingmask), date_predicted)
              featureswbmasked = featureswbmasked.cat(ft)
              
              // RB MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(roadmask).updateMask(buildingmask), date_predicted)
              featuresrbmasked = featuresrbmasked.cat(ft)
              
              // ALL MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(roadmask).updateMask(buildingmask), date_predicted)
              featuresallmasked = featuresallmasked.cat(ft)
            })
        })
        
        var metrics = ee.FeatureCollection(features)
        Export.table.toDrive({
            collection: metrics,
            description: 'STARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featureswatermasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'WSTARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featuresroadsmasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'RSTARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featuresbuildingsmasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'BSTARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featureswrmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '1STARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featureswbmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '2STARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featuresrbmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '3STARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        metrics = ee.FeatureCollection(featuresallmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '4STARFM1_' + roi_abre + '_' + app.picker.pairDate.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
    }



    ////////// Export ////////////////////

    app.exportOriginals = function() {
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())


        if (predictedDate) {

            Export.image.toDrive({
                image: predictedDate.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'L8_Original-' + app.picker.predictedDate.select.getValue(),
                scale: 30,
                region: roi
            });

            var modisPredicted = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', predictedDate.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m']).map(maskCloudsModis)
                .map(addNDVI).first()).resample('bilinear').reproject({
                crs: predictedDate.select('NDVI').projection()
            })

            Export.image.toDrive({
                image: modisPredicted.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'Modis_Original-' + app.picker.predictedDate.select.getValue(),
                scale: 30,
                region: roi
            });
        }
    }

    app.exportPredictedImage = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())

        if (predictedDate && pairDate) {
            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var image = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                // Export the image to Drive.
                Export.image.toDrive({
                    image: image.updateMask(image.select('NDVI').neq(0)).select('NDVI').clip(roi).visualize(ndviParams),
                    description: 'PredictedImage-' + app.picker.predictedDate.select.getValue() + '_W' + app.W_OPTIONS[app.KW.W.select.getValue()].w + 'K' + app.K_OPTIONS[app.KW.K.select.getValue()].k + app.KW.applyNDVI.getValue(),
                    scale: 30,
                    region: roi
                });
            } else {
                // If an image id is found, create an image.
                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate, predictedDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');

                Export.image.toDrive({
                    image: ee.Image(ndvi).updateMask(ndvi.select('NDVI').neq(0)).select('NDVI').clip(roi).visualize(ndviParams),
                    description: 'PredictedImage_NDVIAfter-' + app.picker.predictedDate.select.getValue() + '_W' + app.W_OPTIONS[app.KW.W.select.getValue()].w + 'K' + app.K_OPTIONS[app.KW.K.select.getValue()].k + app.KW.applyNDVI.getValue(),
                    scale: 30,
                    region: roi
                });

            }
        }

    }

    app.exportPairs = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var pairDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())

        if (pairDate) {

            // Export the image to Drive.
            Export.image.toDrive({
                image: pairDate.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'L8_Pair-' + app.picker.pairDate.select.getValue(),
                scale: 30,
                region: roi
            });

            var modisPair = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m']).map(maskCloudsModis)
                .map(addNDVI).first()).resample('bilinear').reproject({
                crs: pairDate.select('NDVI').projection()
            })

            Export.image.toDrive({
                image: modisPair.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'Modis_Pair-' + app.picker.pairDate.select.getValue(),
                scale: 30,
                region: roi
            });
        }


    }


}

/** Creates the app constants. */
app.createConstants = function() {
    app.COLLECTION_ID = 'LANDSAT/LC08/C01/T1_SR';
    app.SECTION_STYLE = {
        margin: '20px 0 0 0'
    };
    app.HELPER_TEXT_STYLE = {
        margin: '8px 0 -3px 8px',
        fontSize: '12px',
        color: 'gray'
    };
    app.IMAGE_COUNT_LIMIT = 25;
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
            nameR: 'Santa Cruz da Trapa e São Cristóvão de Lafões',
            abre: 'SCTSCL',
            roi: saocristovao
        },
        'Benfica do Ribatejo': {
            nameR: 'Benfica do Ribatejo',
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
            description: 'Using 4 Clusters',
            k: 4
        },
        '6': {
            description: 'Using 6 Clusters',
            k: 6
        },
        '8': {
            description: 'Using 8 Clusters',
            k: 8
        }
    }
    app.W_OPTIONS = {
        '3': {
            description: 'Using 3*3 Neighboors',
            w: 3
        },
        '5': {
            description: 'Using 5*5 Neighboors',
            w: 5
        },
        '7': {
            description: 'Using 7*7 Neighboors',
            w: 7
        },
        '9': {
            description: 'Using 9*9 Neighboors',
            w: 9
        }
    }
};

/** Creates the application interface. */
app.boot = function() {
    app.createConstants();
    app.createHelpers();
    app.createPanels();
    var main = ui.Panel({
        widgets: [
            app.intro.panel,
            app.roi.panel,
            app.filters.panel,
            app.picker.panel,
            app.KW.panel,
            app.show.panel,
            app.export.panel,
        ],
        style: {
            width: '570px',
            padding: '8px'
        }
    });
    Map.setCenter(-7.799, 39.20, 6);
    ui.root.insert(1, main);
};

app.boot();