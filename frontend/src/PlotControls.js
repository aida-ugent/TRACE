import { useEffect, useState } from "react";
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
import { DodgedBarplot, StackedBarplot } from "./barplot";
import { explainCluster } from "./api";
import { Explanation } from "./explanation";

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
    exclus,
    selectedPoints,
    children } = props

  const [visibility, setVisibility] = useState('visible')
  const [unstablePointFraction, setUnstablePointFraction] = useState(0.1);
  const [opacityByDensity, setOpacityByDensity] = useState(true);
  const [opacity, setOpacity] = useState({ 'slider': 0.2, 'value': scaleOpacity(0.2) });

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
  }, [visibility])

  const toggleOpacityByDensity = (byDensity) => {
    setOpacityByDensity(byDensity);

    if (scatterplot !== null) {
      if (byDensity) {
        if (scatterplot.get('opacityBy') !== 'valueW') {
          scatterplot.set({
            "opacityBy": "density",
          })
        }
      } else {
        if (scatterplot.get('opacityBy') === 'valueW') {
          scatterplot.set({
            opacityBy: 'w',
            opacity: [opacity['value'], 1],
          })
        } else {
          scatterplot.set({
            opacityBy: null,
            opacity: opacity['value'],
          })
        }
      }
    }
  }

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

  const handleOpacitySelect = (newOpacity) => {
    newOpacity = Math.max(0, Math.min(newOpacity, 1));

    const scaledOpacity = scaleOpacity(newOpacity);
    setOpacity({ 'slider': newOpacity, 'value': scaledOpacity });

    if (scatterplot.get('opacityBy') == 'valueW') {
      scatterplot.set({
        opacityBy: 'w',
        opacity: [scaledOpacity, 1],
      })
    } else {
      scatterplot.set({
        opacityBy: null,
        opacity: scaledOpacity,
      })
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
                {children}

                {/* Point Colors */}
                <div className="flex flex-col items-left my-2 justify-between">
                  <label className="text-sm text-gray-500 w-fit min-w-fit" >point colors</label>
                  <GroupedSelect onChange={pointColorOnChange} options={pointColorOptions} selected={selectedPointColor} />
                </div>

                <div className='flex flex-wrap items-center justify-between my-2'>
                  {/* Point Size */}
                  <div className='flex flex-col w-1/2 items-left pr-2 justify-between'>
                    <label className="text-sm text-gray-500 w-fit min-w-fit" htmlFor="pointSizeSlider">point size</label>
                    <input
                      className="transparent h-[2px] cursor-pointer appearance-none border-transparent bg-neutral-300 mb-2 mt-3"
                      type="range"
                      min={0.1}
                      max={10}
                      step={0.1}
                      value={pointSize}
                      //defaultValue={pointSize}
                      onChange={(event) => handlePointSizeSelect(+event.target.value)}
                      id="pointSizeSlider" />
                  </div>

                  {/* Point Opacity */}
                  <div className='flex flex-col w-1/2 items-left pl-2 justify-between'>
                    <div className="flex flex-row items-center justify-left">
                      <label className="text-sm text-gray-500 w-fit min-w-fit" htmlFor="opacityCheckbox">opacity by density</label>
                      <Checkbox
                        text=""
                        id='opacityCheckbox'
                        checked={opacityByDensity}
                        onChange={toggleOpacityByDensity}
                      />
                    </div>
                    <input
                      className={`transparent h-[2px] w-full appearance-none border-transparent bg-neutral-300 
                  mb-2 mt-3 ${opacityByDensity ? 'accent-slate-100' : 'cursor-pointer'}`}
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

                {/* Infobox */}
                <Infobox
                  scatterplot={scatterplot}
                  selectedPointColor={selectedPointColor}
                  pointColors={pointColors}
                  colorMap={colorMap}
                  pointColorOptions={pointColorOptions} />

                {/* Legend */}
                <div className="flex flex-wrap items-start my-2 justify-start">
                  <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='hoverSwitch'>show legend</label>
                  {/* <div className="w-1/2 flex items-start justify-start"> */}
                  <span className="ml-1">
                    <Switch
                      id='hoverSwitch'
                      checked={legendVisibility == "visible" ? true : false}
                      onChange={toggleLegendVisibility}
                      className={`${legendVisibility == "visible" ? 'bg-blue-600' : 'bg-gray-200'
                        } relative inline-flex h-5 w-9 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                    >
                      <span
                        className={`${legendVisibility == "visible" ? 'translate-x-5' : 'translate-x-1'
                          } inline-block h-3 w-3 transform rounded-full bg-white transition-transform`}
                      />
                    </Switch>
                    {/* </div> */}
                  </span>
                </div>

                {/* Histogram */}
                {pointColors["type"] === "continuous" && pointColors["values"].length > 0 &&
                  <Histogram
                    featureValues={pointColors["values"]}
                    xlabel={selectedPointColor}
                    selectedPoints={selectedPoints}
                    selectedGroupName="selected" />
                }

                {/* Dodged Barplot */}
                {pointColors["type"] === "categorical" &&
                  pointColors["values"].length > 0 &&
                  colorMap["colors"].length > 1 &&
                  <DodgedBarplot
                    featureValues={pointColors["values"]}
                    xlabel={selectedPointColor}
                    selectedPoints={selectedPoints}
                    selectedGroupName="selected"
                    colorMap={colorMap} />
                }

              </Tab>
              <Tab label="Embedding Quality">
                {/* Distance measures */}
                <div className="flex flex-col items-left my-2 justify-between">
                  <div className='flex flex-row'>
                    <label className="text-sm text-gray-500 w-fit min-w-fit mr-1">HD metric</label>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor"
                      className="size-5 cursor-pointer text-gray-500"
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

                <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-3" >
                  High-dimensional neighbors
                </h4>
                <p className="text-sm text-gray-500 text-left my-1">
                  Visualize the high-dimensional neighbors of any point in the 2D embedding to explore the local quality.
                </p>

                {/* K Neighbors */}
                <div className='flex flex-wrap items-start justify-between my-2'>
                  <div className='flex flex-col w-1/2 items-left pr-2 justify-between'>
                    <div className='flex flex-row'>
                      <label className="text-sm text-gray-500 w-fit min-w-fit mr-1" htmlFor="neighborsSlider">neighbors {kNeighbors}</label>
                      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor"
                        className="size-5 cursor-pointer text-gray-500"
                        data-tooltip-id="neighbors-tooltip">
                        <path strokeLinecap="round" strokeLinejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
                      </svg>
                    </div>

                    <input
                      className="transparent h-[2px] cursor-pointer appearance-none border-transparent bg-neutral-300 mb-2 mt-3"
                      type="range"
                      min={0}
                      max={maxNeighbors}
                      step={10}
                      defaultValue={kNeighbors}
                      onChange={(event) => handlekNeighborSelect(+event.target.value)}
                      id="neighborsSlider" />
                  </div>


                  {/* Hover neighbors */}
                  <div className="flex-row items-start text-left justify-left w-1/2 pl-2">
                    <label className="text-sm text-gray-500 w-fit min-w-fit" >show on hover</label>
                    <span className="ml-3">
                      <Switch
                        id='hoverSwitch'
                        checked={hoverNeighborsEnabled}
                        onChange={(enabled) => setHoverNeighborsEnabled(enabled)}
                        className={`${hoverNeighborsEnabled ? 'bg-blue-600' : 'bg-gray-200'
                          } relative inline-flex h-5 w-9 items-center rounded-full transition-colors 
                                    focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                      >
                        <span
                          className={`${hoverNeighborsEnabled ? 'translate-x-5' : 'translate-x-1'
                            } inline-block h-3 w-3 transform rounded-full bg-white transition-transform`}
                        />
                      </Switch>
                    </span>
                  </div>
                </div>

                {/* Compute Neighbors */}
                {/* <AsyncButton onClick={() => precomputeNeighborsSubscribe(maxNeighbors, selectedMetric)}>precompute neighbors</AsyncButton> */}

                <div className='flex flex-wrap items-center mb-2 justify-left'>
                  <span data-tooltip-id='hdneighbors-tooltip'>
                    <AsyncButton onClick={() => handleHDNeighbors(kNeighbors, selectedMetric)}>HD neighbors</AsyncButton>
                  </span>
                  <span data-tooltip-id='intrusions-tooltip'>
                    <AsyncButton onClick={() => showIntrusions(scatterplot, kNeighbors, selectedMetric)}>intrusions</AsyncButton>
                  </span>
                </div>


                <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-3" >
                  High-dimensional distances
                </h4>
                <p className="text-sm text-gray-500 text-left my-1">
                  Select a single point to color points according to their HD distance. The point colors are based on the distances between&nbsp;
                  <a onClick={() => showLandmarks(scatterplot)} className="underline cursor-pointer">landmark points</a>.
                </p>
                <div className='flex flex-wrap items-center mb-2 justify-left'>
                  <DefaultButton onClick={() => pointColorOnChange("HD distances")}>
                    HD distances
                  </DefaultButton>
                </div>

                <span className="select-none text-sm text-gray-500 text-left my-2">
                  <p> Add current point selection to user_annotations.json</p>
                </span>

                <SavePointForm scatterplot={scatterplot} />
              </Tab>

              <Tab label="Exploration">
                <Explanation
                  exclus={exclus}
                  selectedPoints={selectedPoints}
                  pointColorOnChange={pointColorOnChange}
                  pointColorOptions={pointColorOptions} 
                  scatterplot={scatterplot}
                  selectedPointColor={selectedPointColor}/>
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