import { useState } from "react";
import { ReactSelect } from "./utils";
import { SettingsButton, ChevronRightButton, AsyncButton, DefaultButton } from "./buttons";
import { Switch } from '@headlessui/react'
import { Infobox } from "./infobox";
import GroupedSelect from "./groupedSelect";
import { Tabs, Tab } from "./tabs";
import Checkbox, { Radio } from "./checkbox";

const showLandmarks = (scatterplot) => {
  fetch(`/backend/landmarkPoints/`)
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
      fetch("/backend/intrusions", {
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
    fetch(`/backend/precomputeAllNeighbors?maxK=${kNeighbors}&hd_metric=${metric}`, {
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
  let response = await fetch(`/backend/precomputeAllNeighbors?maxK=${maxK}&hd_metric=${metric}`,
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
    children } = props

  const [visibility, setVisibility] = useState('visible')
  const [unstablePointFraction, setUnstablePointFraction] = useState(0.1);
  const [opacityByDensity, setOpacityByDensity] = useState(true);
  const [opacity, setOpacity] = useState({ 'slider': 0.2, 'value': scaleOpacity(0.2) });

  const toggleVisibility = () => {
    if (visibility == "visible") setVisibility("hidden"); else setVisibility("visible");
  }

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
        <div className="select-none right-0 w-1/4 min-w-[300px] h-screen max-h-screen bg-gray-100">
          <span className="relative">
            <div className="absolute -left-[30px] top-2">
              <ChevronRightButton onClick={toggleVisibility} />
            </div>
          </span>

          <div className="flex flex-col h-screen max-h-screen items-top justify-start text-center p-5 pt-12 overflow-y-auto">

            <h3 className="text-lg font-medium leading-6 text-gray-900 w-fit" >
              Settings
            </h3>

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
                  className="transparent h-[2px] cursor-pointer appearance-none 
                border-transparent bg-neutral-300 mb-2 mt-3"
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
              <div className='flex flex-col w-1/2 items-left pr-2 justify-between'>
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
              pointColorOptions={pointColorOptions} />

            {/* Legend */}
            <Switch.Group>
              <div className="flex flex-wrap items-center my-2 justify-between">
                <Switch.Label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='hoverSwitch'>show legend</Switch.Label>
                <div className="w-1/2 flex items-start justify-start">
                  <Switch
                    id='hoverSwitch'
                    checked={legendVisibility == "visible" ? true : false}
                    onChange={toggleLegendVisibility}
                    className={`${legendVisibility == "visible" ? 'bg-blue-600' : 'bg-gray-200'
                      } relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                  >
                    <span
                      className={`${legendVisibility == "visible" ? 'translate-x-6' : 'translate-x-1'
                        } inline-block h-4 w-4 transform rounded-full bg-white transition-transform`}
                    />
                  </Switch>
                </div>
              </div>
            </Switch.Group>


            {/* Neighbors */}
            <h3 className="text-lg font-medium leading-6 text-gray-900 w-fit mt-8" >
              Embedding Quality
            </h3>

            <div className='flex flex-wrap items-start justify-between my-2'>
              <div className='flex flex-col w-1/2 items-left pr-2 justify-between'>
                <label className="text-sm text-gray-500 w-fit min-w-fit" htmlFor="neighborsSlider">neighbors {kNeighbors}</label>
                <input
                  className="transparent h-[2px] cursor-pointer appearance-none 
                border-transparent bg-neutral-300 mb-2 mt-3"
                  type="range"
                  min={0}
                  max={maxNeighbors}
                  step={10}
                  defaultValue={kNeighbors}
                  onChange={(event) => handlekNeighborSelect(+event.target.value)}
                  id="neighborsSlider" />
              </div>

              {/* Distance measures */}
              <div className="flex flex-col w-1/2 items-left pl-2 justify-between">
                <label className="text-sm text-gray-500 w-fit min-w-fit">HD metric</label>
                <ReactSelect
                  options={metricOptions}
                  selected={selectedMetric}
                  onChange={metricOnChange}
                  menuPlacement={'top'}
                />
              </div>
              {/* kNN */}
              {/* <div className='flex flex-col items-left my-2 justify-between'>
              <NamedSlider onChange={handlekNeighborSelect} defaultValue={kNeighbors} min="0" max={maxNeighbors} step="10" label="neighbors" />
            </div> */}
            </div>

            <span className="text-sm text-gray-500 text-left my-2">
              <p>
                Select one or more points</p>
            </span>

            {/* Compute Neighbors */}
            {/* <AsyncButton onClick={() => precomputeNeighborsSubscribe(maxNeighbors, selectedMetric)}>precompute neighbors</AsyncButton> */}

            <div className='flex flex-wrap items-center mb-2 justify-left'>
              <AsyncButton onClick={() => handleHDNeighbors(kNeighbors, selectedMetric)}>HD neighbors</AsyncButton>
              <AsyncButton onClick={() => showIntrusions(scatterplot, kNeighbors, selectedMetric)}>intrusions</AsyncButton>
            </div>

            <span className="select-none text-sm text-gray-500 text-left my-2">
              <p> Select a single point to color points according to their HD distance. The point colors are based on the distances between on&nbsp;
                <a onClick={() => showLandmarks(scatterplot)} className="underline cursor-pointer">landmark points</a>.</p>
            </span>
            <div className='flex flex-wrap items-center mb-2 justify-left'>
              <DefaultButton onClick={() => pointColorOnChange("HD distances")}>
                HD distances
              </DefaultButton>
            </div>

            {/* Unstable points */}
            <h4 className="text-md font-medium leading-6 text-gray-900 w-fit mt-4 flex" >
              Unstable Points
              <abbr title="Highlighting points that change location and neighborhods between the selected embedding and the next.">
                <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="gray" className="w-5 h-5 mx-1">
                  <path fillRule="evenodd" d="M18 10a8 8 0 1 1-16 0 8 8 0 0 1 16 0ZM8.94 6.94a.75.75 0 1 1-1.061-1.061 3 3 0 1 1 2.871 5.026v.345a.75.75 0 0 1-1.5 0v-.5c0-.72.57-1.172 1.081-1.287A1.5 1.5 0 1 0 8.94 6.94ZM10 15a1 1 0 1 0 0-2 1 1 0 0 0 0 2Z" clipRule="evenodd" />
                </svg>
              </abbr>
            </h4>
            <div className="flex flex-wrap items-center my-2 justify-between">
              <label className="text-sm text-gray-500 w-fit min-w-fit" htmlFor="fractionInput">fraction in %</label>
              <input
                className="rounded-md bg-white m-1 px-2 pb-2 pt-2.5 text-sm font-medium leading-normal text-gray-700 "
                type="number"
                min={0.0}
                max={10}
                step={0.001}
                defaultValue={unstablePointFraction}
                onChange={(event) => handleFractionInput(+event.target.value)}
                id="fractionInput" />
              <AsyncButton onClick={() => showUnstablePoints(unstablePointFraction)}>compute</AsyncButton>
            </div>

          </div>
        </div >

        {/*         <Tabs>
          <Tab label="Tab 1">
            <div className="py-4">
              <h2 className="text-lg font-medium mb-2">Tab 1 Content</h2>
              <p className="text-gray-700">Hi Text</p>
            </div>
          </Tab>
          <Tab label="Tab 2">
            <div className="py-4">
              <h2 className="text-lg font-medium mb-2">Tab 2 Content</h2>
              <p className="text-gray-700">Hi Text</p>
            </div>
          </Tab>
        </Tabs> */}
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