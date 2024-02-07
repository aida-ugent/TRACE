import { useEffect, useState } from "react";
import { HoverNote } from "./utils";
import { Switch } from '@headlessui/react'
import GroupedSelect from "./groupedSelect";

function HoverText({ pointId, pointData }) {
    return (
        <>
            {Object.keys(pointData).map((key) => {
                return (
                    <p key={key}>{key}: {pointData[key]}</p>
                )
            })}
        </>
    )
}

export function Infobox(props) {
    const {
        scatterplot,
        selectedPointColor,
        pointColors,
        pointColorOptions,
        ...other } = props

    const [pointOverUnsubscriber, setPointOverUnsubscriber] = useState(null)
    const [pointOutUnsubscriber, setPointOutUnsubscriber] = useState(null)
    const [hoverState, setHoverState] = useState({
        pointId: 0,
        visibility: "hidden",
        position: [0, 0]
    });

    // HOVER FEATURES
    const [hoverNoteEnabled, setHoverNoteEnabled] = useState(false)
    const [hoverData, setHoverData] = useState([])
    const [hoverFeatures, setHoverFeatures] = useState([])
    const hoverFeatureOptions = pointColorOptions.filter((v) => v.label != "quality")
    hoverFeatureOptions.forEach((v) => {
        if (v.label == "metadata") {
            v.options = v.options.filter((v) => v.label != "none" && v.label != "HD distances");
        }
    })

    // when changing the dataset and the pointColorOptions change, we need to reset the hover features
    useEffect(() => {
        setHoverFeatures([]);
        handleHoverNoteEnabled(false)
    }, [pointColorOptions])

    const handlePointOver = (pointId) => {
        setHoverState({
            pointId: pointId,
            visibility: "visible",
            position: scatterplot.getScreenPosition(pointId)
        })
    }

    const hoverFeaturesOnChange = (newValue) => {
        let newFeatures = newValue.map(v => v.value);
        setHoverFeatures(newFeatures);
        if (hoverNoteEnabled) {
            // do we need to fetch new features?
            let fetchFeatures = newFeatures.filter(x => !hoverFeatures.includes(x));
            let deleteFeatures = hoverFeatures.filter(x => !newFeatures.includes(x));
            if (fetchFeatures.length > 0) {
                getHoverData(newFeatures);

                // subscribe to pointover and pointout events
                subscribe();

            } else if (deleteFeatures.length > 0) {
                deleteFeatures.forEach(f => {
                    delete hoverData[f];
                })

                // disable hover not if no features are selected
                if (newFeatures.length == 0) {
                    unsubscribe();
                }
            }
        }
    }

    const prettyPrint = (value) => {
        if (typeof value === 'string') {
            return value;
        } else if (typeof value === 'number') {
            return value.toPrecision(3);
        } else {
            // if it's not a string or a number, just return the value
            return value;
        }
    }

    const handlePointOut = () => {
        setHoverState(prevState => ({
            ...prevState,
            visibility: "hidden"
        }));
    }

    const subscribe = () => {
        // check if we're all subscribed
        if (pointOverUnsubscriber === null) {
            setPointOverUnsubscriber(scatterplot.subscribe('pointover', handlePointOver));
            setPointOutUnsubscriber(scatterplot.subscribe('pointout', handlePointOut));
        }
    }

    const unsubscribe = () => {
        if (pointOverUnsubscriber !== null) {
            scatterplot.unsubscribe(pointOverUnsubscriber)
            scatterplot.unsubscribe(pointOutUnsubscriber)
            setPointOutUnsubscriber(null);
            setPointOverUnsubscriber(null);
        }
    }

    const getHoverColor = (pointId) => {
        if (pointColors["type"] === "categorical") {
            return pointColors["colorMap"][pointColors["values"][pointId]];
        } else {
            return ("#e9e9e9");
        }
    }

    const getHoverData = (feature) => {
        let body = JSON.stringify({ feature_list: feature })
        return new Promise((resolve, reject) => {
            if (feature.length == 0) {
                resolve(true);
            } else {
                console.log("fetching hover metadata features")
                fetch("/backend/metadataFeatures", {
                    method: "POST",
                    body: body,
                    headers: {
                        "Content-type": "application/json; charset=UTF-8"
                    }
                })
                    .then(res => res.json())
                    .then(data => {
                        setHoverData(data);
                        resolve(true);
                    })
            }
        })
    }

    const getSingleHoverData = (pointID) => {
        // make a dictionary with keys the hover features and values the hover data for point with index pointID
        if (hoverState.visibility == "visible") {
            let pointData = {}
            hoverFeatures.forEach(f => {
                pointData[f] = prettyPrint(hoverData[f][pointID])
            })
            return pointData;
        } else {
            return {}
        }
    }

    const handleHoverNoteEnabled = (enabled) => {
        if (scatterplot !== null) {
            if (enabled) {
                if (hoverFeatures.length > 0) {
                    getHoverData(hoverFeatures)
                        .then(() => subscribe())
                }
            }
            else unsubscribe();
        }
        setHoverNoteEnabled(enabled);
    }

    return (
        <>
            <HoverNote
                visible={hoverState.visibility}
                position={hoverState.position}
                color={getHoverColor(hoverState.pointId)}>
                <HoverText
                    pointId={hoverState.pointId}
                    pointColors={pointColors}
                    selectedPointColor={selectedPointColor}
                    pointData={getSingleHoverData(hoverState.pointId)} />
            </HoverNote>
            <div className="flex flex-col items-left my-2 justify-between">
                <Switch.Group>
                    <div className="flex items-start justify-start">
                        <label className="text-sm text-gray-500 w-fit min-w-fit" >show features on hover</label>
                        <span className="ml-3 mb-2">
                            <Switch
                                id='hoverSwitch'
                                checked={hoverNoteEnabled}
                                onChange={handleHoverNoteEnabled}
                                className={`${hoverNoteEnabled ? 'bg-blue-600' : 'bg-gray-200'
                                    } relative inline-flex h-5 w-9 items-center rounded-full transition-colors 
                                    focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                            >
                                <span
                                    className={`${hoverNoteEnabled ? 'translate-x-5' : 'translate-x-1'
                                        } inline-block h-3 w-3 transform rounded-full bg-white transition-transform`}
                                />
                            </Switch>
                        </span>
                    </div>
                </Switch.Group>
                <GroupedSelect
                    onChange={hoverFeaturesOnChange}
                    options={hoverFeatureOptions}
                    selected={hoverFeatures}
                    isMulti={true}
                    isDisabled={!hoverNoteEnabled}
                    isClearable={true} />
            </div>
        </>
    )
}