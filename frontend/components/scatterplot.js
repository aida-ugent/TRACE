'use client'

import { useCallback, useState, useEffect } from "react";
import { DefaultButton, ResetButton, AsyncButton } from "./buttons";
import Legend from "./legend";
import { scatterplot } from "./canvas";
import { SettingsMenu } from "./PlotControls"
import SelectionInfo from "./selectionInfo";
import { Colorbar } from "./colorbar";
import CanvasWrapper from '@/components/canvas_wrapper'
import { ReactSelect } from "./utils"
import GroupedSelect from "./groupedSelect"
import { getHDNeighbors } from "./api";

// detault options
let pointSizeInitial = 5;
let maxNeighbors = 200;

let filteredPoints = [];
let numPoints;

const resetZoomHandler = () => {
    scatterplot.zoomToLocation([0.0, 0.0], 1.2,
        {
            transition: true,
            transitionDuration: 1000,
        });
}

const resetOpacityHandler = () => {
    scatterplot.set({
        opacityBy: 'density',
    })
}

const zoomToSelectionHandler = (opacities) => {
    const selected = scatterplot.get('selectedPoints')
    opacities.forEach((v, i) => { if (v == 1) { selected.push(i) } })
    console.log(`zooming to ${selected.length} points`)
    scatterplot.zoomToPoints(selected, {
        padding: 1.0,
        transition: true,
        transitionDuration: 1500,
    })
}

const flattenGroupedOptions = (groupedOptions) => {
    return groupedOptions.map(obj => obj.options.map(v => v.value)).flat();
}

const filterPoints = (featureValues, filterValue) => {
    // only show points with value != 'filterValue'
    let show = [];
    filteredPoints.map((v) => {
        if (featureValues[v] != filterValue) {
            show.push(v)
        }
    });
    let hiddenLength = filteredPoints.length - show.length;
    filteredPoints = show;
    scatterplot.filter(show);
    console.log(`hiding ${filterValue} (${hiddenLength} points)`)
}

const unfilterPoints = (featureValues, filterValue) => {
    // show points with value 'filterValue'
    let show = [];
    featureValues.map((v, i) => {
        if (v == filterValue) show.push(i);
    })
    let addedLength = show.length;
    show = show.concat(filteredPoints)
    filteredPoints = show;
    scatterplot.filter(show)
    console.log(`showing ${filterValue} (${addedLength} points)`)
}

const resetPointFilter = () => {
    filteredPoints = [...Array(numPoints).keys()];
}

const getPointColors = (embName, featureName) => {
    var fetchStr = `/backend/pointColor/${featureName}?embeddingName=${embName}`;

    if (scatterplot != null) {
        const selected = scatterplot.get('selectedPoints')
        if (selected.length == 1) {
            fetchStr = `/backend/pointColor/${featureName}?embeddingName=${embName}&selectedPoint=${selected[0]}`;
        }
    }

    return new Promise((resolve, reject) => {
        fetch(fetchStr)
            .then(res => res.json())
            .then(res => {
                resolve(res)
            })
    })
}


export const isSameElements = (a, b) => {
    if (a.length !== b.length) return false;
    const aSet = new Set(a);
    return b.every((value) => aSet.has(value));
};

const showEmbedding = (embedding, pointColor, opacities, useTransition = true, preventFilterReset = true, opacityBy = null) => {
    let cview = scatterplot.get('cameraView');
    if (!preventFilterReset) resetPointFilter();

    scatterplot
        .draw(
            {
                x: embedding['x'],
                y: embedding['y'],
                z: pointColor["encoded_values"],
                w: opacities
            },
            {
                transition: useTransition,
                transitionDuration: 1500,
                transitionEasing: 'quadInOut',
                // preventFilterReset: preventFilterReset,
                filter: preventFilterReset ? filteredPoints : undefined,
                zDataType: pointColor["type"]
            }
        ).then(() => {
            scatterplot.set({
                cameraView: cview,
                pointColor: Object.values(pointColor["colorMap"]),
                opacityBy: opacityBy != null ? opacityBy : scatterplot.get('opacityBy'),
                opacity: [0.03, 1],
                //pointColorActive: "#ffba08", //"#55308d"
            })
        })
}


const fetchEmbedding = (embName) => {
    return new Promise((resolve, reject) => {
        fetch(`/backend/embedding?embName=${embName}`)
            .then(res => res.json())
            .then(res => {
                resolve(res)
            })
    })
}


export default function Scatterplot() {
    const [isLoading, setIsLoading] = useState(true);
    const [loadingMessage, setLoadingMessage] = useState("");
    const [isLoadingData, setIsLoadingData] = useState(false);
    const [isError, setIsError] = useState(false);
    const [errorMessage, setErrorMessage] = useState("");

    const [datasetOptions, setDatasetOptions] = useState([]);
    const [embeddingName, setEmbeddingName] = useState("UMAP");
    const [activeEmbedding, setActiveEmbedding] = useState(null)
    const [embeddingOptions, setEmbeddingOptions] = useState(null);

    const [legendVisibility, setLegendVisibility] = useState("visible");
    const [kNeighbors, setkNeighbors] = useState(50);
    const [pointColorOptions, setPointColorOptions] = useState(null);
    const [selectedPointColor, setSelectedPointColor] = useState("initialPointColors");
    const [pointColors, setPointColors] = useState(
        {
            "values": 0,
            "encoded_values": 0,
            "colorMap": { "none": "#444444" },
            "type": "categorical",
            "group": "metadata"
        }
    );
    const [opacities, setOpacities] = useState(null);
    const [scatterLoaded, setScatterLoaded] = useState(false);
    const [scatterplotState, setScatterplot] = useState(null);
    const [datasetInfo, setDatasetInfo] = useState("");
    const [selectedMetric, setMetric] = useState("angular");
    const [metricOptions, setMetricOptions] = useState([]);
    const [datasetName, setDatasetName] = useState(null);

    const handleDatasetSelect = (newDatasetName) => {
        console.log(`handleDatasetSelect ${newDatasetName}`)
        setDatasetName(newDatasetName);
        loadDatasetConfiguration(newDatasetName, setIsLoadingData);
    }

    const loadDatasetConfiguration = (datasetName, setLoadingFunction) => {
        setLoadingFunction(true);
        return new Promise((resolve, reject) => {
            var newEmbeddingName;
            setLoadingMessage(`Loading ${datasetName}...`)
            fetch(`/backend/loadDataset?datasetName=${datasetName}`)
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/loadDataset?datasetName=${datasetName} `
                            + `got HTTP error: Status ${response.status}, ${response.statusText}`);
                        setIsError(true);
                        setErrorMessage(`Error loading dataset ${datasetName}. Refresh page to start again.`)
                        reject(response)
                    } else {
                        response.json()
                            .then(res => {
                                setMetric(res["hd_metric"]);
                                setMetricOptions(res["metric_options"].map(v => { return { value: v, label: v } }));
                                setDatasetInfo(res["dataset_info"]);
                                setDatasetName(res["dataset_name"]);
                                setEmbeddingOptions(res["embedding_options"])
                                let newPointColorOptions = [];
                                let groupOptions = Object.keys(res["point_color_options"])
                                let newSelectedPointColor = res["point_color_options"][groupOptions[0]][0]
                                for (const key of groupOptions) {
                                    newPointColorOptions.push({
                                        'label': key,
                                        'options': res["point_color_options"][key].map(v => { return { value: v, label: v } })
                                    })

                                }
                                setPointColorOptions(newPointColorOptions);
                                setSelectedPointColor(newSelectedPointColor);
                                console.log(newPointColorOptions)

                                newEmbeddingName = flattenGroupedOptions(res["embedding_options"])[0];
                                setEmbeddingName(newEmbeddingName);

                                setLoadingMessage(`Fetching ${newEmbeddingName} embedding...`)
                                Promise.all([fetchEmbedding(newEmbeddingName), getPointColors(newEmbeddingName, newSelectedPointColor)])
                                    .then(([newEmbedding, pointColors]) => {
                                        setActiveEmbedding(newEmbedding)
                                        numPoints = newEmbedding["x"].length;
                                        filteredPoints = [...Array(numPoints).keys()];
                                        setOpacities(new Array(numPoints).fill(0));
                                        console.log(`number of points ${numPoints}`)
                                        setPointColors(pointColors);
                                        console.log("Finished loading data")
                                        setLoadingFunction(false);
                                        setLoadingMessage("");
                                        resolve(true);
                                    })
                            })
                    }
                })
                .catch(err => {
                    console.log(err)
                    setIsError(true);
                    setErrorMessage(`Error loading dataset ${datasetName}.`)
                    reject(err)
                }
                );
        })
    }

    const getNextEmbeddingName = (currEmbeddingName) => {
        let embeddingNames = flattenGroupedOptions(embeddingOptions);
        let currIndex = embeddingNames.indexOf(currEmbeddingName);
        let nextIndex = (currIndex + 1) % embeddingNames.length;
        return embeddingNames[nextIndex];
    }

    const getPreviousEmbeddingName = (currEmbeddingName) => {
        let embeddingNames = flattenGroupedOptions(embeddingOptions);
        let currIndex = embeddingNames.indexOf(currEmbeddingName);
        let nextIndex = (currIndex - 1 + embeddingNames.length) % embeddingNames.length;
        return embeddingNames[nextIndex];
    }

    const handlePointColorSelect = (newPointColor) => {
        console.log(`handlePointColorSelect ${newPointColor}`)
        setSelectedPointColor(newPointColor);
        getPointColors(embeddingName, newPointColor)
            .then((res) => {
                setPointColors(res);
                showEmbedding(activeEmbedding, res, opacities, false, false);
            })
    }

    const handleHDNeighbors = (kNeighbors, metric) => {
        let selectedPoints = scatterplot.get('selectedPoints');
        return new Promise((resolve, reject) => {
            if (selectedPoints.length > 0) {
                getHDNeighbors(selectedPoints, kNeighbors, metric)
                    .then((binary_neighbors) => {
                        setOpacities(binary_neighbors);
                        showEmbedding(activeEmbedding, pointColors, binary_neighbors, false, true, "w");
                        resolve(true);
                    })
            } else {
                resolve(true);
            }
        })
    }

    // fetch the list of unstable points from the backend
    const showUnstablePoints = (maxFraction) => {
        return new Promise((resolve, reject) => {
            scatterplot.deselect();
            fetch(`/backend/getUnstablePoints?embNameA=${embeddingName}` +
                `&embNameB=${getNextEmbeddingName(embeddingName)}` +
                `&maxFraction=${maxFraction}` +
                `&k=${kNeighbors}`)
                .then(res => res.json())
                .then(res => {
                    console.log(`showing ${res['result'].length} unstable points`)
                    scatterplot.select(res['result']);
                    resolve(true);
                })
        })
    }

    const handleEmbeddingSelect = (newEmbeddingName) => {
        setEmbeddingName(newEmbeddingName);
        fetchEmbedding(newEmbeddingName)
            .then(newEmbedding => {
                setActiveEmbedding(newEmbedding);

                if (pointColors["group"] === "quality") {
                    // pointColor was quality of old embedding ... recompute
                    getPointColors(newEmbeddingName, selectedPointColor)
                        .then((res) => {
                            if ("none" in res["colorMap"]) {
                                setSelectedPointColor("none");
                            }
                            setPointColors(res);
                            showEmbedding(newEmbedding, res, opacities);
                        })
                } else {
                    showEmbedding(newEmbedding, pointColors, opacities);
                }
            })
    }

    useEffect(() => {
        fetch("/backend/datasetOptions")
            .then(response => {
                if (!response.ok) {
                    console.log(`/backend/datasetOptions got HTTP error: ` +
                        `Status ${response.status}, Text ${response.statusText}`);
                    setIsError(true);
                    setErrorMessage(`Error loading dataset options.`)
                } else {
                    response.json()
                        .then(res => {
                            setDatasetName(res["result"][0])
                            setDatasetOptions(res["result"]);
                            loadDatasetConfiguration(res["result"][0], setIsLoadingData)
                                .then(res => setIsLoading(false), err => setIsLoading(true))
                        })
                }
            })
    }, []);

    useEffect(() => {
        if (scatterLoaded && !isLoading && !isLoadingData) {
            console.log(`drawing initial embedding ${embeddingName}`)
            showEmbedding(activeEmbedding, pointColors, opacities, false, false);
            scatterplot.deselect();
            resetOpacityHandler();
        }
    }, [scatterLoaded, isLoading, isLoadingData]);

    if (isError) {
        return (
            <div role="status" className="absolute w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center flex-col">
                <svg xmlns="http://www.w3.org/2000/svg"
                    fill="none"
                    viewBox="0 0 24 24"
                    strokeWidth={1.5}
                    stroke="currentColor"
                    className="w-10 h-10">
                    <path strokeLinecap="round" strokeLinejoin="round" d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" />
                </svg>
                <p className="text-gray-500 text-sm">{errorMessage}</p>
            </div>
        );
    } else if (isLoading & !isError) {
        return (
            <div role="status" className="absolute w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center flex-col">
                <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600 m-4" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                <p className="text-gray-500 text-sm">{loadingMessage}</p>
            </div>
        );
    } else {
        return (
            <>
                <div className="flex-grow h-screen max-h-screen bg-white relative min-w-[400px]">
                    <CanvasWrapper setScatterLoaded={setScatterLoaded} setScatterplot={setScatterplot} />
                    <div className="absolute w-full flex flex-wrap top-0 pt-1 px-1 items-center display-block justify-between">
                        <ResetButton onClick={resetZoomHandler} />

                        {/* Embedding method */}
                        <div className="flex items-center p-2 justify-between">
                            <button
                                type="button"
                                className="h-fit select-none inline-flex rounded-full bg-slate-500 mx-2 p-1 text-sm font-medium 
                                leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 
                                focus:shadow-none focus:outline-none hover:bg-slate-800"
                                onClick={() => handleEmbeddingSelect(getPreviousEmbeddingName(embeddingName))}
                            >
                                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={3} stroke="currentColor" className="w-5 h-5">
                                    <path strokeLinecap="round" strokeLinejoin="round" d="M15.75 19.5 8.25 12l7.5-7.5" />
                                </svg>
                            </button>
                            <GroupedSelect onChange={handleEmbeddingSelect} options={embeddingOptions} selected={embeddingName} />
                            <button
                                type="button"
                                className="h-fit select-none inline-flex rounded-full bg-slate-500 mx-2 p-1 text-sm font-medium 
                                leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 
                                focus:shadow-none focus:outline-none hover:bg-slate-800"
                                onClick={() => handleEmbeddingSelect(getNextEmbeddingName(embeddingName))}
                            >
                                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={3} stroke="currentColor" className="w-5 h-5">
                                    <path strokeLinecap="round" strokeLinejoin="round" d="m8.25 4.5 7.5 7.5-7.5 7.5" />
                                </svg>
                            </button>
                        </div>
                        <div style={{ visibility: 'hidden' }} />
                    </div>
                    <SelectionInfo scatterplot={scatterplotState} numPoints={numPoints} datasetInfo={datasetInfo} />
                </div>
                <SettingsMenu
                    scatterplot={scatterplotState}
                    selectedPointColor={selectedPointColor}
                    pointColors={pointColors}
                    pointSizeInitial={pointSizeInitial}
                    pointColorOptions={pointColorOptions}
                    pointColorOnChange={handlePointColorSelect}
                    setLegendVisibility={setLegendVisibility}
                    legendVisibility={legendVisibility}
                    kNeighbors={kNeighbors}
                    maxNeighbors={maxNeighbors}
                    handlekNeighborSelect={setkNeighbors}
                    metricOptions={metricOptions}
                    selectedMetric={selectedMetric}
                    metricOnChange={setMetric}
                    showUnstablePoints={showUnstablePoints}
                    handleHDNeighbors={handleHDNeighbors}
                >
                    {/* Dataset */}
                    <div className="flex flex-col items-left my-2 justify-between">
                        <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" >Data</label>
                        <ReactSelect
                            options={datasetOptions.map(v => { return { value: v, label: v } })}
                            selected={datasetName}
                            onChange={handleDatasetSelect}
                        />
                    </div>

                </SettingsMenu>
                <div className="fixed bottom-0 left-0 flex-wrap flex-row m-2 max-h-[90%] overflow-auto">
                    {
                        pointColors["type"] === "categorical" ?
                            <Legend
                                colormap={pointColors["colorMap"]}
                                title={selectedPointColor}
                                filter={(filterValue) => filterPoints(pointColors["values"], filterValue)}
                                unfilter={(filterValue) => unfilterPoints(pointColors["values"], filterValue)}
                                visibility={legendVisibility}
                            /> :
                            <Colorbar
                                colormap={pointColors["colorMap"]}
                                title={selectedPointColor}
                                visibility={legendVisibility}
                            />
                    }
                </div>


                <div
                    //visibility={isLoadingData ? "visible" : "hidden"}
                    role="status" className="absolute w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center flex-col"
                    style={{ visibility: isLoadingData ? "visible" : "hidden" }}
                >
                    <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600 m-4" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                    <span className="text-gray-500 text-sm">{loadingMessage}</span>
                </div>
            </>
        );
    }
}


