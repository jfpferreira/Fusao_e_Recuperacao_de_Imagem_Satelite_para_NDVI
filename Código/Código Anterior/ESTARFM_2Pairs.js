var modis = ee.ImageCollection("MODIS/006/MOD09GQ"),
    L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    osmwater = ee.FeatureCollection("users/jfpferreira/OSM_Water"),
    osmroads = ee.FeatureCollection("users/jfpferreira/OSM_Roads"),
    osmbuildings = ee.FeatureCollection("users/jfpferreira/OSM_Buildings");
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
fine1 - Landsat Image at date 1 (pair)
fine2 - Landsat Image at date 2 (pair)
lst_predicted - original image Landsat, used to compare the predicted image to the original
Window_pixels - Window, can be 3,5,7,9
K - use for clustering. Can be 4,6,8
band - NIR, RED or NDVI
ROI - From the 4 available
*/
var predict = function(window_pixels, k, band, roi, fine1, fine2, lst_predicted) {

    var modis_pair1 = ee.Image(modis
        .filterBounds(roi)
        .map(addTime)
        .filterMetadata('Date', 'equals', fine1.get('Date'))
        .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
        .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
        crs: fine1.select(band).projection()
    })

    var modis_pair2 = ee.Image(modis
        .filterBounds(roi)
        .map(addTime)
        .filterMetadata('Date', 'equals', fine1.get('Date'))
        .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
        .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
        crs: fine2.select(band).projection()
    })


    var modis_predicted = ee.Image(modis
        .filterBounds(roi)
        .map(addTime)
        .filterMetadata('Date', 'equals', lst_predicted.get('Date'))
        .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
        .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
        crs: lst_predicted.select(band).projection()
    })



    var f1 = ee.ImageCollection(fine1).select([band, 'latitude', 'longitude', 'Time'])
        .map(function(img) {
            return img.addBands(modis_pair1.select(band))
                .addBands(modis_predicted.select(band))
                .rename([band + '_LSTPair', 'latitude', 'longitude', 'Time', band + '_ModisPair', band + '_ModisPredicted'])
                .set('Time', img.get('system:time_start'))
        })

    var f2 = ee.ImageCollection(fine2).select([band, 'latitude', 'longitude', 'Time'])
        .map(function(img) {
            return img.addBands(modis_pair2.select(band))
                .addBands(modis_predicted.select(band))
                .rename([band + '_LSTPair', 'latitude', 'longitude', 'Time', band + '_ModisPair', band + '_ModisPredicted'])
                .set('Time', img.get('system:time_start'))
        })

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
        }).get('groups')) //.aside(print,'Range Clusters '+band)

        return img.addBands(result, ['cluster'])

    }

    var candidates = function(name, img) {

        return img.select(name).neighborhoodToBands(
            ee.Kernel.square({
                radius: ee.Number(window_pixels).subtract(1).divide(2),
                units: 'pixels'
            })).toArray().rename(name + '_candidates')
    }

    var similarPixels = function(img) {

        var allCandidatesArrays = img.addBands(candidates('cluster', img))
            .addBands(candidates(band + '_LSTPair', img))
            .addBands(candidates(band + '_ModisPair', img))
            .addBands(candidates(band + '_ModisPredicted', img))
            .addBands(candidates('slope', img))
            .addBands(candidates('latitude', img))
            .addBands(candidates('longitude', img))

        var mask_clusters = (allCandidatesArrays.select('cluster').eq(allCandidatesArrays.select('cluster_candidates'))).and(allCandidatesArrays.select('latitude'))
            .and(allCandidatesArrays.select(band + '_LSTPair').neq(ee.Image(-2)))

        var image = allCandidatesArrays.select(band + '_LSTPair_candidates', band + '_ModisPair_candidates', band + '_ModisPredicted_candidates', 'cluster_candidates', 'slope_candidates', 'latitude_candidates', 'longitude_candidates').arrayMask(mask_clusters)

        var nCandidates = image.select('cluster_candidates').arrayLengths().arrayGet([0]).rename('nCandidates')

        return img.addBands(image).addBands(nCandidates.select('nCandidates'))

    }

    var weight = function(img, original) {

        var correlation = img.select([band + '_LSTPair_candidates', band + '_ModisPair_candidates']).toArray(1)
            .arrayReduce(ee.Reducer.pearsonsCorrelation(), [0], 1)
            .arrayProject([1]).arrayFlatten([
                ['c', 'p']
            ])

        var euclDist = img.expression('sqrt(pow((lt-ltc),2) + pow((lg-lgc),2))', {
            ltc: img.select(['latitude_candidates']),
            lgc: img.select(['longitude_candidates']),
            lt: original.select(['latitude']),
            lg: original.select(['longitude'])
        }).rename('distance')

        var count_candidates = euclDist.select('distance').arrayReduce(ee.Reducer.count(), [0]).arrayGet([0])
        var distance = ee.Image(1).divide(ee.Image(1).add(euclDist.select('distance')).divide(count_candidates.divide(2))).rename('distance')
        var sum_distance = distance.select('distance').arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0])

        return distance.select('distance').divide(sum_distance).rename('weight')

    }

    var conversionCoeff = function(fine, coarse, coarsePredicted) {

        // It is a linear regression
        //Independent: Time
        //Dependent: NDVI

        var fine_ndvi = ee.ImageCollection(fine.select([band + '_LSTPair', 'Time'])).map(function(img) {
            return img.rename([band, 'Time'])
        })

        var sortedCollection = ee.ImageCollection(coarse).merge(ee.ImageCollection(coarsePredicted)).merge(fine_ndvi).filterBounds(roi).sort('Time');

        var linearFit = sortedCollection.select(['Time', band]) //Without damage pixels
            .reduce(ee.Reducer.linearFit())

        return linearFit.select(["scale"]).rename('slope')

    }

    var predictByImage = function(img, weight, conversionCoeff) {

        var coarseSub = img.select(band + '_ModisPredicted_candidates').subtract(img.select(band + '_ModisPair_candidates')).rename('coarseSub')

        var weightConver = weight.select('weight').multiply(conversionCoeff.select('slope_candidates')).multiply(coarseSub.select('coarseSub'))
            .arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0]).rename('weightConver')

        return img.select(band + '_LSTPair').add(weightConver.select('weightConver')).rename(band) //predictedNDVI

    }

    var temporalWeight = function(img) {

        return ee.Image(1).divide(img.select(band + '_ModisPredicted_candidates').subtract(img.select(band + '_ModisPair_candidates'))
            .arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0]).abs())

    }

    // Linear regression to each pixel
    var coef1 = conversionCoeff(f1, modis_pair1, modis_predicted)
    var coef2 = conversionCoeff(f2, modis_pair2, modis_predicted)

    // Clustering including masked pixels and add to them slope band
    var cluster_fine1 = clustering(ee.Image(f1.first()).unmask(-2)).addBands(coef1);
    var cluster_fine2 = clustering(ee.Image(f2.first()).unmask(-2)).addBands(coef2);

    //Add neighbourhood and their data
    var img_fine1 = similarPixels(cluster_fine1)
    var img_fine2 = similarPixels(cluster_fine2)

    //distance and correlation
    var weight1 = weight(img_fine1, cluster_fine1)
    var weight2 = weight(img_fine2, cluster_fine2)

    //Predict pixels for each image
    var predict1 = predictByImage(img_fine1, weight1, img_fine1)
    var predict2 = predictByImage(img_fine2, weight2, img_fine2)

    var final1 = temporalWeight(img_fine1).divide(temporalWeight(img_fine1).add(temporalWeight(img_fine2))).multiply(predict1).rename(band)
    var final2 = temporalWeight(img_fine2).divide(temporalWeight(img_fine1).add(temporalWeight(img_fine2))).multiply(predict2).rename(band)

    var predicted = final1.add(final2).addBands(img_fine1.select('nCandidates').add(img_fine2.select('nCandidates')))
        .addBands(modis_predicted.select(band).rename(band + '_ModisPredicted')).addBands(ee.Image.pixelLonLat())
        .addBands(lst_predicted.select(band).rename(band + '_LSTPredicted'))

    return predicted.addBands(lst_predicted.select([band], [band + '_LSTPredicted']).subtract(predicted.select(band)).rename('error_LST'))
        .addBands(modis_predicted.select([band], [band + '_ModisPredicted']).subtract(predicted.select(band)).rename('error_Modis'))
        .addBands(cluster_fine1.select(['cluster'], ['cluster1'])).addBands(cluster_fine2.select(['cluster'], ['cluster2']))

}

///////////////////////////////// Analysis //////////////////////////////////////////////////////
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

var histograms = function(date_pair1, date_pair2, date_predicted, roi, roi_name, flag) {

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
                var pre = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted)

                var histogram = ui.Chart.image.histogram(pre.updateMask(pre.neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).select(['NDVI', 'NDVI_LSTPredicted']), roi, 30)
                    .setSeriesNames(['NDVI_Predicted', 'NDVI_LST'])
                    .setOptions(options);

                // Display the histogram.
                print(histogram);
            })

            var k4 = predict(w[indexW], 4, 'NDVI', roi, date_pair1, date_pair2, date_predicted).select(['NDVI'], ['NDVI_K4']).updateMask(date_predicted.select('NDVI'))
            var k6 = predict(w[indexW], 6, 'NDVI', roi, date_pair1, date_pair2, date_predicted).select(['NDVI'], ['NDVI_K6']).updateMask(date_predicted.select('NDVI'))
            var k8 = predict(w[indexW], 8, 'NDVI', roi, date_pair1, date_pair2, date_predicted).select(['NDVI'], ['NDVI_K8']).addBands(k4).addBands(k6).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2)).updateMask(k4.neq(0))

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

        var image = predict(3, 4, 'NDVI', roi, date_pair1, date_pair2, date_predicted).updateMask(date_predicted.select('NDVI'))

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

                var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_pair2, date_predicted);
                var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_pair2, date_predicted);
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
            var imageNIR4 = predict(w[indexW], 4, 'nir', roi, date_pair1, date_pair2, date_predicted);
            var imageRED4 = predict(w[indexW], 4, 'red', roi, date_pair1, date_pair2, date_predicted);
            var k4 = imageNIR4.select('nir').addBands(imageRED4.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI_K4').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
            // K=6
            var imageNIR6 = predict(w[indexW], 6, 'nir', roi, date_pair1, date_pair2, date_predicted);
            var imageRED6 = predict(w[indexW], 6, 'red', roi, date_pair1, date_pair2, date_predicted);
            var k6 = imageNIR6.select('nir').addBands(imageRED6.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI_K6').updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))

            // K=8
            var imageNIR8 = predict(w[indexW], 8, 'nir', roi, date_pair1, date_pair2, date_predicted);
            var imageRED8 = predict(w[indexW], 8, 'red', roi, date_pair1, date_pair2, date_predicted);
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

        var image = predict(3, 4, 'NDVI', roi, date_pair1, date_pair2, date_predicted).updateMask(date_predicted.select('NDVI'))

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

var scatters = function(date_pair1, date_pair2, date_predicted, roi, roi_name, flag) {

    var w = [3, 5, 7, 9];
    var k = [4, 6, 8];

    if (flag) {
        w.forEach(function(nameW, indexW) {
            k.forEach(function(nameK, indexK) {

                var predicted = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted)
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

                var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_pair2, date_predicted);
                var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_pair2, date_predicted);
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


/////////////////////////////////////// Visualization ///////////////////////////////////////////

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
            app.KW.K.select,
            app.KW.W.select,
            app.filters.applyButton,
            app.picker.predictedDate.select,
            app.picker.pairDate1.select,
            app.picker.pairDate2.select
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
            app.picker.pairDate1.select.items().reset(ids);
            app.picker.pairDate2.select.items().reset(ids);
            // Default the image picker to the first id.
            app.picker.predictedDate.select.setValue(app.picker.predictedDate.select.items().get(0));
            app.picker.pairDate1.select.setValue(app.picker.pairDate1.select.items().get(0));
            app.picker.pairDate2.select.setValue(app.picker.pairDate2.select.items().get(0));
        });

        collection = filtered;
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
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first());
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first());

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
            maps[0].addLayer(originalDate.select('NDVI').clip(roi), ndviParams, 'OriginalImage');
            maps[0].addLayer(originalDate.select('NDVI').unmask(-2).eq(-2).updateMask(originalDate.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');

        }

        if (pairDate1 && pairDate2 && originalDate) {

            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var predicted = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', roi,
                    pairDate1, pairDate2, originalDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                maps[1].addLayer(predicted.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                maps[1].addLayer(predicted.select('NDVI').eq(0).updateMask(predicted.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)
                print(predicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

            } else {

                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', roi,
                    pairDate1, pairDate2, originalDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', roi,
                    pairDate1, pairDate2, originalDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                maps[1].addLayer(ndvi.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                maps[1].addLayer(ndvi.select('NDVI').eq(0).updateMask(ndvi.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
                maps[1].add(legend)
                print(ndvi.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')

            }

        }
        ui.root.widgets().set(0, mapGrid);

    };

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
            map.addLayer(predictedDate.select('NDVI').unmask(-1).eq(-1).updateMask(predictedDate.select('NDVI').unmask(-1).eq(-1)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(predictedDate.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original')

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

        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first());

        var modisPredicted = ee.Image(modis
            .filterBounds(roi)
            .map(addTime)
            .filterMetadata('Date', 'equals', predictedDate.get('Date'))
            .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
            .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
            crs: predictedDate.select('nir').projection()
        })

        map.addLayer(modisPredicted.select('NDVI').clip(roi), ndviParams, 'OriginalModisPredicted');
        print(modisPredicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Original Modis')
    };

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
        map.add(ui.Label('Predicted Image'));
        map.add(legend);
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        map.centerObject(roi, 12);
        map.addLayer(roi, {}, 'ROI');
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())


        if (pairDate1 && pairDate2 && predictedDate) {
            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var image = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', roi,
                    pairDate1, pairDate2, predictedDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                map.addLayer(image.select('NDVI').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), ndviParams, 'PredictedImage');
                map.addLayer(image.select('NDVI').eq(0).updateMask(image.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');
            } else {
                // If an image id is found, create an image.
                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, predictedDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, predictedDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                map.addLayer(ndvi.select('NDVI').clip(roi), ndviParams, 'PredictedImage');
                map.addLayer(ndvi.select('NDVI').eq(0).updateMask(ndvi.select('NDVI').eq(0)).clip(roi), {
                    palette: '#b30000'
                }, 'masked');

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
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())


        var NAMES = [
            'Pair1 LST',
            'Pair1 Modis',
            'Pair2 LST',
            'Pair2 Modis'
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
        maps[3].setControlVisibility({
            scaleControl: true
        });

        // Create a grid of maps.
        var mapGrid = ui.Panel(
            [
                ui.Panel([maps[0], maps[1]], null, {
                    stretch: 'both'
                }),
                ui.Panel([maps[2], maps[3]], null, {
                    stretch: 'both'
                })
            ],
            ui.Panel.Layout.Flow('horizontal'), {
                stretch: 'both'
            }
        );

        if (pairDate1 && pairDate2) {

            maps[0].addLayer(pairDate1.select('NDVI').clip(roi), ndviParams, 'Pair1LST');
            maps[0].addLayer(pairDate1.select('NDVI').unmask(-2).eq(-2).updateMask(pairDate1.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(pairDate1.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair1 LST')

            var modisPair1 = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate1.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
                .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
                crs: pairDate1.select('nir').projection()
            })

            print(modisPair1.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair1 Modis')
            maps[1].addLayer(modisPair1.select('NDVI').clip(roi), ndviParams, 'PairModis1');
            maps[1].add(legend);

            // Pair2

            maps[2].addLayer(pairDate2.select('NDVI').clip(roi), ndviParams, 'Pair2LST');
            maps[2].addLayer(pairDate2.select('NDVI').unmask(-2).eq(-2).updateMask(pairDate2.select('NDVI').unmask(-2).eq(-2)).clip(roi), {
                palette: '#b30000'
            }, 'masked');
            print(pairDate2.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair2 LST')

            var modisPair2 = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate2.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
                .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
                crs: pairDate2.select('nir').projection()
            })

            print(modisPair2.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Pair2 Modis')
            maps[3].addLayer(modisPair2.select('NDVI').clip(roi), ndviParams, 'PairModis2');
            //maps[3].add(legend);
        }

        ui.root.widgets().set(0, mapGrid);

    };

    app.refreshMapClusters = function() {

        Map.clear();
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var k = app.K_OPTIONS[app.KW.K.select.getValue()].k
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())

        var NAMES = [
            'Clusters LST1',
            'Clusters LST2'
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


        if (app.KW.applyNDVI.getValue()) {
            if (pairDate1 && pairDate2 && predictedDate) {
                // If an image id is found, create an image.
                maps[0].addLayer(clusters(pairDate1, 'NDVI', k, roi).select('cluster').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi).randomVisualizer(), {}, 'Clusters-PairImage1');
                maps[1].addLayer(clusters(pairDate2, 'NDVI', k, roi).select('cluster').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi).randomVisualizer(), {}, 'Clusters-PairImage2');
            }
        } else {
            print('Impossible to Show Clusters. Select apply NDVI from the beginning')
        }
        ui.root.widgets().set(0, mapGrid);
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
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())

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

        if (pairDate1 && pairDate2 && originalDate) {

            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var predicted = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, originalDate);
                // Add the image to the map with the corresponding visualization options.
                //var visOption = app.VIS_OPTIONS[app.vis.select.getValue()];
                maps[0].addLayer(predicted.select('error_Modis').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), errorParams, 'Error_Modis');
                maps[1].addLayer(predicted.select('error_LST').clip(app.ROI_OPTIONS[app.roi.select.getValue()].roi), errorParams, 'Error_LST');
                maps[1].add(legend)
                print(predicted.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')
            } else {

                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', roi,
                    pairDate1, pairDate2, originalDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', roi,
                    pairDate1, pairDate2, originalDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                var error_LST = originalDate.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')


                var modisPredicted = ee.Image(modis
                    .filterBounds(roi)
                    .map(addTime)
                    .filterMetadata('Date', 'equals', originalDate.get('Date'))
                    .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m'])
                    .map(addNDVI).map(maskCloudsModis).first()).resample('bilinear').reproject({
                    crs: originalDate.select('B5').projection()
                })

                var error_Modis = modisPredicted.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_Modis')
                maps[1].addLayer(error_LST.select('error_LST').clip(roi), errorParams, 'Error_LST');
                maps[0].addLayer(error_Modis.select('error_Modis').clip(roi), errorParams, 'Error_Modis');
                maps[1].add(legend)
                print(error_LST.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean error_LST')
                print(error_Modis.reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean error_Modis')

            }
        }
        ui.root.widgets().set(0, mapGrid);
    };

    app.showHistograms = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var date_pair2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())


        print('--------- Histograms--------------')
        histograms(date_pair1, date_pair2, date_predicted, roi, roi_name, app.KW.applyNDVI.getValue());
    }

    app.showScatters = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var date_pair2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())

        print('--------- Scatter-----------------')
        scatters(date_pair1, date_pair2, date_predicted, roi, roi_name, app.KW.applyNDVI.getValue());

    }

    app.metrics = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var date_pair2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())
        app.exportMetrics()
        var w = [3, 5, 7, 9];
        var k = [4, 6, 8];

        if (app.KW.applyNDVI.getValue()) {

            print('----------------R Squared NDVI-----------------')
            w.forEach(function(nameW, indexW) {
                k.forEach(function(nameK, indexK) {

                    print('---------K=' + nameK + '------------')
                    var image = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted)
                    image = image.updateMask(image.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))

                    print(rsquared(date_predicted.updateMask(image.select('NDVI').neq(0)), image, 'NDVI', roi), 'R squared W' + nameW + 'K' + nameK)

                    // Original - Predicted
                    var error = date_predicted.select('NDVI').subtract(image.select('NDVI')).rename('error')
                    var error_dict = ee.Dictionary(['Error', error.reduceRegion(ee.Reducer.mean(), roi, 30).get('error')])
                    var mae = ee.Dictionary(['MAE', error.select('error').abs().reduceRegion(ee.Reducer.mean(), roi, 30).get('error')]);
                    var rmse = ee.Dictionary(['RMSE', ee.Number((error.select('error').multiply(error.select('error'))).select('error').reduceRegion(ee.Reducer.mean(), roi, 30).get('error')).sqrt()])

                    var dict = mae.combine(rmse).combine(error_dict);
                    print(dict, 'Stats W' + nameW + 'K' + nameK)
                    print(predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted).updateMask(image.select('NDVI').neq(0)).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')


                    print(image.select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30), 'pearsons Correlation');
                    print(predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted).updateMask(image.select('NDVI').neq(0)).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted Image W' + nameW + 'K' + nameK)

                })
            })

        } else {
            print('----------- R Squared Reflectance--------------')
            w.forEach(function(nameW, indexW) {
                k.forEach(function(nameK, indexK) {

                    print('---------K=' + nameK + '------------')
                    var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_pair2, date_predicted);
                    var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_pair2, date_predicted);
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

                    print(image.updateMask(image.select('NDVI').neq(0)).reduceRegion(ee.Reducer.mean(), roi, 30), 'Mean Predicted')


                    print(image.select(['NDVI', 'NDVI_LSTPredicted']).reduceRegion(ee.Reducer.pearsonsCorrelation(), roi, 30), 'pearsons Correlation');

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

    ////////// Export ////////////////////
    app.exportOriginals = function() {
        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())

        if (predictedDate) {
            // Export the image to Drive.
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
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())

        if (predictedDate && pairDate1 && pairDate2) {
            if (app.KW.applyNDVI.getValue()) {

                // If an image id is found, create an image.
                var image = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'NDVI', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, predictedDate);
                // Export the image to Drive.
                Export.image.toDrive({
                    image: ee.Image(image).updateMask(image.select('NDVI').neq(0)).select('NDVI').clip(roi).visualize(ndviParams),
                    description: 'PredictedImage-' + app.picker.predictedDate.select.getValue(),
                    scale: 30,
                    region: roi
                });

                // Export the image to Drive.
                Export.image.toDrive({
                    image: ee.Image(image).select('error_LST').clip(roi).visualize(errorParams),
                    description: 'Error-' + app.picker.predictedDate.select.getValue(),
                    scale: 30,
                    region: roi
                });


            } else {
                // If an image id is found, create an image.
                var imageNIR = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'nir', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, predictedDate);
                var imageRED = predict(app.W_OPTIONS[app.KW.W.select.getValue()].w, app.K_OPTIONS[app.KW.K.select.getValue()].k, 'red', app.ROI_OPTIONS[app.roi.select.getValue()].roi,
                    pairDate1, pairDate2, predictedDate);
                var ndvi = imageNIR.select('nir').addBands(imageRED.select('red')).normalizedDifference(['nir', 'red']).rename('NDVI');
                var error_LST = predictedDate.select('NDVI').subtract(ndvi.select('NDVI')).rename('error_LST')

                Export.image.toDrive({
                    image: ee.Image(ndvi).updateMask(ndvi.select('NDVI').neq(0)).select('NDVI').clip(roi).visualize(ndviParams),
                    description: 'PredictedImage_NDVIAfter-' + app.picker.predictedDate.select.getValue(),
                    scale: 30,
                    region: roi
                });

                // Export the image to Drive.
                Export.image.toDrive({
                    image: ee.Image(error_LST).select('error_LST').clip(roi).visualize(errorParams),
                    description: 'Error-' + app.picker.predictedDate.select.getValue(),
                    scale: 30,
                    region: roi
                });

            }
        }

    }

    app.exportPairs = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var predictedDate = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var pairDate1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var pairDate2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())

        if (pairDate1 && pairDate2) {
            // Export the image to Drive.
            Export.image.toDrive({
                image: pairDate1.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'L8_Pair1-' + app.picker.pairDate1.select.getValue(),
                scale: 30,
                region: roi
            });

            Export.image.toDrive({
                image: pairDate2.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'L8_Pair2-' + app.picker.pairDate2.select.getValue(),
                scale: 30,
                region: roi
            });

            // Modis
            var modisPair1 = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate1.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m']).map(maskCloudsModis)
                .map(addNDVI).first()).resample('bilinear').reproject({
                crs: pairDate1.select('NDVI').projection()
            })

            Export.image.toDrive({
                image: modisPair1.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'Modis_Pair1-' + app.picker.pairDate1.select.getValue(),
                scale: 30,
                region: roi
            });

            var modisPair2 = ee.Image(modis
                .filterBounds(roi)
                .map(addTime)
                .filterMetadata('Date', 'equals', pairDate2.get('Date'))
                .select(['sur_refl_b02', 'sur_refl_b01', 'Time', 'QC_250m'], ['nir', 'red', 'Time', 'QC_250m']).map(maskCloudsModis)
                .map(addNDVI).first()).resample('bilinear').reproject({
                crs: pairDate2.select('NDVI').projection()
            })

            Export.image.toDrive({
                image: modisPair2.select('NDVI').clip(roi).visualize(ndviParams),
                description: 'Modis_Pair2-' + app.picker.pairDate2.select.getValue(),
                scale: 30,
                region: roi
            });
        }


    }

    app.exportMetrics = function() {

        var roi = app.ROI_OPTIONS[app.roi.select.getValue()].roi
        var roi_abre = app.ROI_OPTIONS[app.roi.select.getValue()].abre
        var roi_name = app.ROI_OPTIONS[app.roi.select.getValue()].nameR
        var date_predicted = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.predictedDate.select.getValue()).first())
        var date_pair1 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
        var date_pair2 = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())

        var w = [3, 5, 7, 9];
        var k = [4, 6, 8];
        var dici = ee.Dictionary();
        var features = ee.List([]);
        var featuresbuildingsmasked = ee.List([]);
        var featureswatermasked = ee.List([]);
        var featuresroadsmasked = ee.List([]);
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
        
                
        var makePrediction = function(indexW, indexK, nameW, nameK, date_pair1, date_pair2, date_predicted) {
          var feats = ee.List([]);
          // NDVI
          var imageNDVI = predict(w[indexW], k[indexK], 'NDVI', roi, date_pair1, date_pair2, date_predicted)
          
          var image = imageNDVI.updateMask(imageNDVI.select('NDVI').neq(0)).updateMask(date_predicted.select('NDVI').unmask(-2).neq(-2))
          
          var ft = calc_Metrics(date_predicted, image, nameW, nameK)
          feats = feats.cat(ee.List([ft]))
          
          //return ft;

          // Reflectance
          var imageNIR = predict(w[indexW], k[indexK], 'nir', roi, date_pair1, date_pair2, date_predicted);
          var imageRED = predict(w[indexW], k[indexK], 'red', roi, date_pair1, date_pair2, date_predicted);         
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
              var ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1, date_pair2, date_predicted)
              features = features.cat(ft)
              
              // WATER MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask), date_pair2.updateMask(watermask), date_predicted)
              featureswatermasked = featureswatermasked.cat(ft)
              
              // ROAD MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(roadmask), date_pair2.updateMask(roadmask), date_predicted)
              featuresroadsmasked = featuresroadsmasked.cat(ft)
              
              // BUILDING MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(buildingmask), date_pair2.updateMask(buildingmask), date_predicted)
              featuresbuildingsmasked = featuresbuildingsmasked.cat(ft)
              
              // WR MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(roadmask), date_pair2.updateMask(watermask).updateMask(roadmask), date_predicted)
              featureswrmasked = featureswrmasked.cat(ft)
              
              // WB MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(buildingmask), date_pair2.updateMask(watermask).updateMask(buildingmask), date_predicted)
              featureswbmasked = featureswbmasked.cat(ft)
              
              // RB MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(roadmask).updateMask(buildingmask), date_pair2.updateMask(roadmask).updateMask(buildingmask), date_predicted)
              featuresrbmasked = featuresrbmasked.cat(ft)
              
              // ALL MASK
              ft = makePrediction(indexW, indexK, nameW, nameK, date_pair1.updateMask(watermask).updateMask(roadmask).updateMask(buildingmask), date_pair2.updateMask(watermask).updateMask(roadmask).updateMask(buildingmask), date_predicted)
              featuresallmasked = featuresallmasked.cat(ft)
            })
        })
        
        var metrics = ee.FeatureCollection(features)
        Export.table.toDrive({
            collection: metrics,
            description: 'ESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        var metrics = ee.FeatureCollection(featureswatermasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'WESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        var metrics = ee.FeatureCollection(featuresroadsmasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'RESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        var metrics = ee.FeatureCollection(featuresbuildingsmasked)
        Export.table.toDrive({
            collection: metrics,
            description: 'BESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        var metrics = ee.FeatureCollection(featureswrmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '1ESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
                
        var metrics = ee.FeatureCollection(featureswbmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '2ESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
        
        var metrics = ee.FeatureCollection(featuresrbmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '3ESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
                
        var metrics = ee.FeatureCollection(featuresallmasked)
        Export.table.toDrive({
            collection: metrics,
            description: '4ESTARFM2_' + roi_abre + '_' + app.picker.pairDate1.select.getValue() + '_' +app.picker.pairDate2.select.getValue() + '_' + app.picker.predictedDate.select.getValue(),
            fileFormat: 'CSV',
            folder: 'metrics2',
            selectors: ['W', 'K', 'signal', 'R^2', 'RMSE', 'MAE', 'Error', 'correlation']
        })
    }
};
/** Creates the app constants. */
/** Creates the UI panels. */
app.createPanels = function() {

    app.intro = {
        panel: ui.Panel([
            ui.Label({
                value: 'ESTARFM 2Pairs',
                style: {
                    fontWeight: 'bold',
                    fontSize: '24px',
                    margin: '10px 5px',
                    shown: true
                }
            }),
            ui.Label('This app allows you to apply ESTARFM using 2Pairs Landsat-Modis')
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
        pairDate1: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate1.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.pairDate1.label.setValue('Cloud Cover in ROI: ' + result.toFixed(2))
                    })
                }
            })
        },
        pairDate2: {
            label: ui.Label(),
            select: ui.Select({
                placeholder: 'Select Image ID',
                onChange: function() {
                    var image = ee.Image(collection.filterMetadata('system:index', 'equals', app.picker.pairDate2.select.getValue()).first())
                    var cloud = image.get('CLOUD_COVER_REGION')
                    cloud.evaluate(function(result) {
                        app.picker.pairDate2.label.setValue('Cloud Cover in ROI: ' + result.toFixed(2))
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
                ui.Label('Pair Date1', app.HELPER_TEXT_STYLE), app.picker.pairDate1.select,
                app.picker.pairDate1.label
            ], ui.Panel.Layout.flow('horizontal')),
            ui.Panel([
                ui.Label('Pair Date2', app.HELPER_TEXT_STYLE), app.picker.pairDate2.select,
                app.picker.pairDate2.label
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
        buttonOriginalModis: ui.Button({
            label: 'Show Original Modis',
            // React to the button's click event.
            onClick: function() {
                app.refreshMapOriginalModis()
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

    /* The panel for the export section with corresponding widgets. */
    app.show.panel = ui.Panel({
        widgets: [
            ui.Label('4) Visualization', {
                fontWeight: 'bold'
            }),
            ui.Panel([
                ui.Label('Compare', app.HELPER_TEXT_STYLE),
                app.show.buttonPairs,
                app.show.buttonOriginalPredicted
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
                app.show.buttonScatters
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
            app.export.buttonPredicted,
            app.export.buttonOriginals,
            app.export.buttonPairs
        ],
        style: app.SECTION_STYLE
    });

}

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
            nameR: 'Santa da Cruz da Trapa e São Cristóvão de Lafões',
            abre: 'SCTSCL',
            roi: saocristovao
        },
        'Benfica do Ribatejo ': {
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
            width: '458px',
            padding: '8px'
        }
    });
    Map.setCenter(-7.799, 39.20, 6);
    ui.root.insert(1, main);
};

app.boot();