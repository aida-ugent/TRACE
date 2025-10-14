import { useEffect, useState, useMemo, useCallback } from "react";
import { ReactSelect, SavePointForm } from "./utils";
import { SettingsButton, ChevronRightButton, AsyncButton, DefaultButton } from "./buttons";
import { Switch } from '@headlessui/react'
import { Infobox } from "./infobox";
import GroupedSelect from "./groupedSelect";
import { Tabs, Tab } from "./tabs";
import Checkbox, { Radio } from "./checkbox";
import { backend_url } from "./api";
import { Tooltip } from 'react-tooltip';
import { Histogram } from "./histogram";
import { RangeFilter } from "./RangeFilter";
import { DodgedBarplot, StackedBarplot } from "./barplot";

const showLandmarks = (scatterplot) => {
  fetch(`${backend_url}/backend/landmarkPoints/`)
    .then(response => {
      if (response.ok) {
        response.json()
          .then(res => {
            scatterplot.deselect();
            scatterplot.select(res["result"]);
          })
      }
    })
}


const showIntrusions = (scatterplot, kNeighbors, metric) => {
  let selectedPoints = scatterplot.get('selectedPoints')

  return new Promise((resolve, reject) => {
    if (selectedPoints.length > 1) {
      fetch(`${backend_url}/backend/intrusions`, {
        method: "POST",
        body: JSON.stringify({
          k: kNeighbors,
          points: selectedPoints,
          hd_metric: metric
        }),
        headers: {
          "Content-type": "application/json; charset=UTF-8"
        }
      })
        .then(res => res.json())
        .then(data => {
          let intrusions = data["result"]
          scatterplot.deselect()
          if (intrusions.length > 0) {
            scatterplot.select(intrusions);
            scatterplot.zoomToPoints(selectedPoints, {
              padding: 0.9,
              transition: true,
              transitionDuration: 1500,
            })
          }
          console.log(`showing ${intrusions.length} intrusions`)
          resolve(true);
        })
    } else {
      resolve(true);
    }
  })
}

const precomputeNeighbors = (kNeighbors, metric) => {
  return new Promise((resolve, reject) => {
    fetch(`${backend_url}/backend/precomputeAllNeighbors?maxK=${kNeighbors}&hd_metric=${metric}`, {
      method: "POST",
      body: {},
      headers: {
        "Content-type": "application/json; charset=UTF-8"
      }
    })
      .then(response => {
        if (!response.ok) {
          if (response.status == 500) {
            console.log(`HTTP error: Status ${response.status}, Text ${response.statusText}`);
            return precomputeNeighbors(maxK, metric);
          } else {
            reject(response);
          }

        } else {
          resolve(true);
        }
      })
  })
}

async function precomputeNeighborsSubscribe(maxK, metric) {
  let response = await fetch(`${backend_url}/backend/precomputeAllNeighbors?maxK=${maxK}&hd_metric=${metric}`,
    {
      method: "POST",
      body: {},
      headers: {
        "Content-type": "application/json; charset=UTF-8"
      }
    });

  if (!response.ok) {
    console.log(`HTTP error: Status ${response.status}, Text ${response.statusText}`);
    await new Promise(resolve => setTimeout(resolve, 1000));
    await precomputeNeighborsSubscribe(maxK, metric);
  } else {
    console.log("Precomputed neighbors successfully.");
  }
}


function scaleOpacity(sliderValue) {
  return Math.pow(sliderValue, 3);
}

export function SettingsMenu(props) {
  const {
    scatterplot,
    selectedPointColor,
    pointColors,
    colorMap,
    pointColorOptions,
    pointColorOnChange,
    pointSize,
    setPointSize,
    setLegendVisibility,
    legendVisibility,
    kNeighbors,
    maxNeighbors,
    handlekNeighborSelect,
    metricOptions,
    selectedMetric,
    metricOnChange,
    showUnstablePoints,
    handleHDNeighbors,
    pointColorScaling,
    handlePointColorScaling,
    hoverNeighborsEnabled,
    setHoverNeighborsEnabled,
    selectedPoints,
    opacityByDensity,
    toggleOpacityByDensity,
    opacity,
    handleOpacitySelect,
    rangeFilterResetTrigger,
    handleRangeFilterChange,
    children } = props

  const [visibility, setVisibility] = useState('visible')
  const [unstablePointFraction, setUnstablePointFraction] = useState(0.1);
  
  // Range filter state - persists across tab switches
  // Current filter values (slider handles)
  const [rangeFilterMin, setRangeFilterMin] = useState(null);
  const [rangeFilterMax, setRangeFilterMax] = useState(null);

  // Maximal range (min/max of all values)
  const [maximalRangeMin, setMaximalRangeMin] = useState(null);
  const [maximalRangeMax, setMaximalRangeMax] = useState(null);

  // Initialize range filter to full range when pointColors changes
  useEffect(() => {
    if (pointColors && pointColors["type"] === "continuous" && pointColors.values && pointColors.values.length > 0) {
      // Calculate maximal range
      const numericValues = pointColors.values.filter(val => typeof val === 'number' && !isNaN(val));
      if (numericValues.length > 0) {
        const min = Math.min(...numericValues);
        const max = Math.max(...numericValues);
        if (isFinite(min) && isFinite(max) && min !== max) {
          setMaximalRangeMin(min);
          setMaximalRangeMax(max);
          setRangeFilterMin(min);
          setRangeFilterMax(max);
        }
      }
    } else {
      // Reset for non-continuous data
      setMaximalRangeMin(null);
      setMaximalRangeMax(null);
      setRangeFilterMin(null);
      setRangeFilterMax(null);
    }
    handleRangeFilterChange(null, null);
  }, [pointColors]);

  const toggleVisibility = () => {
    if (visibility == "visible") setVisibility("hidden"); else setVisibility("visible");
  }

  useEffect(() => {
    const canvas = document.getElementById("canvas");
    if (canvas !== null) {
      const { width, height } = canvas.getBoundingClientRect();
      console.log(`canvas width: ${width}, height: ${height}`)
      //scatterplot.set({ width, height });
    }
  }, [visibility]);

  // Local range filter handler that stores values and calls parent
  const handleRangeChange = (minVal, maxVal) => {
    // Only handle range changes for continuous data
    if (!pointColors || pointColors["type"] !== "continuous") {
      return;
    }
    setRangeFilterMin(minVal);
    setRangeFilterMax(maxVal);
    // Compare to maximal range
    const isMinFiltered = (minVal !== null && maximalRangeMin !== null) ? minVal > maximalRangeMin : false;
    const isMaxFiltered = (maxVal !== null && maximalRangeMax !== null) ? maxVal < maximalRangeMax : false;
    if (isMinFiltered || isMaxFiltered) {
      handleRangeFilterChange(isMinFiltered ? minVal : null, isMaxFiltered ? maxVal : null);
    } else {
      handleRangeFilterChange(null, null);
    }
  };

  const handlePointSizeSelect = (pointSize) => {
    setPointSize(pointSize);
    scatterplot.set({ pointSize });
  };

  const handleFractionInput = (fraction) => {
    setUnstablePointFraction(fraction);
  }

  const toggleLegendVisibility = () => {
    if (legendVisibility == "visible") {
      setLegendVisibility("hidden")
    } else {
      setLegendVisibility("visible")
    }
  }


  if (visibility == "visible") {
    return (
      <>
        <div className="select-none right-0 w-1/4 min-w-[380px] h-screen max-h-screen bg-gray-100">
          <div className="relative">
            <div className="absolute -left-[32px]">
              <ChevronRightButton onClick={toggleVisibility} />
            </div>
          </div>

          <div className="flex flex-col h-screen max-h-screen items-top justify-start text-center overflow-y-auto">
            <Tabs>
              <Tab label="Settings">
                <div className="px-4 py-2 space-y-4">
                  {children}

                  {/* Point Colors */}
                  <div className="space-y-2 text-left">
                    <label className="text-sm text-gray-600 font-medium">Point Colors</label>
                    <GroupedSelect 
                      onChange={pointColorOnChange} 
                      options={pointColorOptions} 
                      selected={selectedPointColor} 
                    />
                  </div>

                  <div className="grid grid-cols-2 gap-4 text-left">
                    {/* Point Size */}
                    <div className="space-y-2">
                      <label className="text-sm text-gray-600 font-medium" htmlFor="pointSizeSlider">
                        Point Size
                      </label>
                      <div className="px-1">
                        <input
                          className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer slider"
                          type="range"
                          min={0.1}
                          max={10}
                          step={0.1}
                          value={pointSize}
                          onChange={(event) => handlePointSizeSelect(+event.target.value)}
                          id="pointSizeSlider" 
                        />
                      </div>
                    </div>

                    {/* Point Opacity */}
                    <div className="space-y-2">
                      <div className="flex items-center justify-between">
                        <label className="text-sm text-gray-600 font-medium" htmlFor="opacityCheckbox">
                          Opacity
                        </label>
                        <div className="flex items-center space-x-2">
                          <span className="text-xs text-gray-500">Auto</span>
                          <Checkbox
                            text=""
                            id='opacityCheckbox'
                            checked={opacityByDensity}
                            onChange={toggleOpacityByDensity}
                          />
                        </div>
                      </div>
                      <div className="px-1">
                        <input
                          className={`w-full h-2 bg-gray-200 rounded-lg appearance-none slider ${
                            opacityByDensity ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer'
                          }`}
                          type="range"
                          min={0}
                          max={1}
                          step={0.001}
                          value={opacity['slider']}
                          onChange={(event) => handleOpacitySelect(+event.target.value)}
                          id="pointOpacitySlider"
                          disabled={opacityByDensity}
                        />
                      </div>
                    </div>
                  </div>

                  {/* Legend Toggle */}
                  <div className="space-y-2 text-left">
                    <div className="flex items-center justify-between py-2 px-3 bg-gray-50 rounded-lg">
                      <label className="text-sm text-gray-600" htmlFor='legendSwitch'>
                        show legend
                      </label>
                      <Switch
                        id='legendSwitch'
                        checked={legendVisibility == "visible" ? true : false}
                        onChange={toggleLegendVisibility}
                        className={`${legendVisibility == "visible" ? 'bg-blue-600' : 'bg-gray-300'
                          } relative inline-flex h-5 w-9 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2`}
                      >
                        <span
                          className={`${legendVisibility == "visible" ? 'translate-x-5' : 'translate-x-1'
                            } inline-block h-3 w-3 transform rounded-full bg-white transition-transform shadow-sm`}
                        />
                      </Switch>
                    </div>
                  </div>

                  {/* Infobox */}
                  <div className="space-y-2 text-left">
                    <div className="bg-gray-50 rounded-lg p-3">
                      <Infobox
                        scatterplot={scatterplot}
                        selectedPointColor={selectedPointColor}
                        pointColors={pointColors}
                        colorMap={colorMap}
                        pointColorOptions={pointColorOptions} 
                      />
                    </div>
                  </div>

                  {/* Histogram */}
                  {pointColors["type"] === "continuous" && pointColors["values"].length > 0 &&
                    <div className="space-y-2 text-left">
                      <Histogram
                        featureValues={pointColors["values"]}
                        xlabel={selectedPointColor}
                        selectedPoints={selectedPoints}
                        selectedGroupName="selected" 
                      />
                      
                      {/* Range Filter */}
                      <div className="pt-2">
                        <RangeFilter
                          featureValues={pointColors["values"]}
                          title={selectedPointColor}
                          onRangeChange={handleRangeChange}
                          minValue={rangeFilterMin}
                          maxValue={rangeFilterMax}
                          fullRangeMin={maximalRangeMin}
                          fullRangeMax={maximalRangeMax}
                        />
                      </div>
                    </div>
                  }

                  {/* Dodged Barplot */}
                  {pointColors["type"] === "categorical" &&
                    pointColors["values"].length > 0 &&
                    colorMap["colors"].length > 1 &&
                    <div className="space-y-2 text-left">
                      <DodgedBarplot
                        featureValues={pointColors["values"]}
                        xlabel={selectedPointColor}
                        selectedPoints={selectedPoints}
                        selectedGroupName="selected"
                        colorMap={colorMap} 
                      />
                    </div>
                  }
                </div>

              </Tab>
              <Tab label="Embedding Quality">
                <div className="px-4 py-2 space-y-4">
                  {/* HD Metric Selection */}
                  <div className="space-y-2 text-left">
                    <div className='flex flex-row items-center'>
                      <label className="text-sm text-gray-600 font-medium mr-1">HD Metric</label>
                      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor"
                        className="size-4 text-gray-400 cursor-pointer"
                        data-tooltip-id="hdmetric-tooltip">
                        <path strokeLinecap="round" strokeLinejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
                      </svg>
                    </div>
                    <ReactSelect
                      options={metricOptions}
                      selected={selectedMetric}
                      onChange={metricOnChange}
                      menuPlacement={'top'}
                    />
                  </div>

                  {/* High-dimensional neighbors section */}
                  <div className="space-y-2 text-left">
                    <h4 className="text-sm text-gray-600 font-medium">
                      High-dimensional neighbors
                    </h4>
                    <p className="text-sm text-gray-500 mb-3">
                      Visualize the high-dimensional neighbors of any point in the 2D embedding to explore the local quality.
                    </p>

                    <div className="grid grid-cols-2 gap-4">
                      {/* K Neighbors */}
                      <div className="space-y-2">
                        <div className='flex flex-row items-center'>
                          <label className="text-sm text-gray-600 font-medium mr-1" htmlFor="neighborsSlider">
                            Neighbors ({kNeighbors})
                          </label>
                          <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor"
                            className="size-4 text-gray-400 cursor-pointer"
                            data-tooltip-id="neighbors-tooltip">
                            <path strokeLinecap="round" strokeLinejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
                          </svg>
                        </div>
                        <div className="px-1">
                          <input
                            className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer slider"
                            type="range"
                            min={0}
                            max={maxNeighbors}
                            step={10}
                            value={kNeighbors}
                            onChange={(event) => handlekNeighborSelect(+event.target.value)}
                            id="neighborsSlider" 
                          />
                        </div>
                      </div>

                      {/* Hover neighbors */}
                      <div className="space-y-2">
                        <label className="text-sm text-gray-600 font-medium">Show on hover</label>
                        <div className="px-1 -mt-1">
                          <Switch
                            id='hoverSwitch'
                            checked={hoverNeighborsEnabled}
                            onChange={(enabled) => setHoverNeighborsEnabled(enabled)}
                            className={`${hoverNeighborsEnabled ? 'bg-blue-600' : 'bg-gray-300'
                              } relative inline-flex h-5 w-9 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2`}
                          >
                            <span
                              className={`${hoverNeighborsEnabled ? 'translate-x-5' : 'translate-x-1'
                                } inline-block h-3 w-3 transform rounded-full bg-white transition-transform shadow-sm`}
                            />
                          </Switch>
                        </div>
                      </div>
                    </div>

                    {/* Action Buttons */}
                    <div className='flex flex-wrap gap-2 pt-2'>
                      <span data-tooltip-id='hdneighbors-tooltip'>
                        <AsyncButton onClick={() => handleHDNeighbors(kNeighbors, selectedMetric)}>HD neighbors</AsyncButton>
                      </span>
                      <span data-tooltip-id='intrusions-tooltip'>
                        <AsyncButton onClick={() => showIntrusions(scatterplot, kNeighbors, selectedMetric)}>intrusions</AsyncButton>
                      </span>
                    </div>
                  </div>

                  {/* High-dimensional distances section */}
                  <div className="space-y-2 text-left">
                    <h4 className="text-sm text-gray-600 font-medium">
                      High-dimensional distances
                    </h4>
                    <p className="text-sm text-gray-500 mb-3">
                      Select a single point to color points according to their HD distance. The point colors are based on the distances between&nbsp;
                      <a onClick={() => showLandmarks(scatterplot)} className="underline cursor-pointer text-blue-600 hover:text-blue-700">landmark points</a>.
                    </p>
                    <div>
                      <DefaultButton onClick={() => pointColorOnChange("HD distances")}>
                        HD distances
                      </DefaultButton>
                    </div>
                  </div>

                  {/* Save annotations section */}
                  <div className="space-y-2 text-left">
                    <h4 className="text-sm text-gray-600 font-medium">Save annotations</h4>
                    <p className="text-sm text-gray-500 mb-3">
                      Add current point selection to user_annotations.json
                    </p>
                    <SavePointForm scatterplot={scatterplot} />
                  </div>
                </div>
              </Tab>

            </Tabs>

          </div>
        </div >

        {/* TOOLTIPS */}
        <Tooltip id="neighbors-tooltip" className='max-w-[300px] text-sm text-left'>
          Choose how many of the HD neighbors (based on {selectedMetric} distance) will be shown.
          For this dataset a maximum of {maxNeighbors} can be selected.
        </Tooltip>
        <Tooltip id="hdmetric-tooltip" className='max-w-[300px] text-sm text-left'>
          Select the metric that was used to precompute the HD neighbors and distances between points.
        </Tooltip>
        <Tooltip id="hdneighbors-tooltip" className='max-w-[300px] text-sm text-left'>
          Select a group of points to show the union of their HD neighbors. Use the lasso while holding shift or select several points with Ctrl.
        </Tooltip>
        <Tooltip id="intrusions-tooltip" className='max-w-[300px] text-sm text-left'>
          Select a group of points to show intrusions. These are points which are selected but are not within the union of HD neighbors of the selected points.
        </Tooltip>
      </>
    )
  } else {
    return (
      <>
        <div className="fixed right-0 top-0 my-2">
          <SettingsButton onClick={toggleVisibility} />
        </div>
      </>
    )
  }
}