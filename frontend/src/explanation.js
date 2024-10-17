import { useEffect, useState } from "react";
import { explainCluster, compareClusters } from "./api";
import { DefaultButton } from "./buttons";
import { ReactSelect } from "./utils";


const selectPointsScatterplot = (selectedPoints, scatterplot) => {
    scatterplot.deselect();
    scatterplot.select(selectedPoints);
}


export const ComparingClusters = ({ selectedPoints, pointColorOnChange, pointColorOptions, scatterplot }) => {
    const [explanation, setExplanation] = useState({
        "features": [],
        "selectionA": [],
        "selectionB": []
    });

    const [selectionA, setSelectionA] = useState([]);
    const [selectionB, setSelectionB] = useState([]);


    const saveSelectionA = (selectedPoints) => {
        setSelectionA(selectedPoints);
    };

    const saveSelectionB = (selectedPoints) => {
        setSelectionB(selectedPoints);
    };

    useEffect(() => {
        setExplanation({
            "features": [],
            "selectionA": [],
            "selectionB": []
        });
    }, [pointColorOptions]);

    const compareClusterExplanation = (selectionA, selectionB) => {
        if (selectionA.length <= 1 || selectionB.length <= 1) {
            return;
        }
        compareClusters(selectionA, selectionB).then((explanation) => {
            const newExplanation = {
                "features": explanation,
                "selectionA": selectionA,
                "selectionB": selectionB
            };
            setExplanation(newExplanation);
        });
    }

    return (
        <span>
            <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-3" >
                Cluster vs. Cluster
            </h4>
            <p className="text-sm text-gray-500 text-left my-1">
                Select two groups of points to compute features that stand out in one group compared to the other.
                If points are in both groups, they will be removed from selection A.
            </p>
            <div className='flex flex-wrap items-center mb-2 justify-start'>
                <DefaultButton onClick={() => saveSelectionA(selectedPoints)}>
                    update selection A ({selectionA.length})
                </DefaultButton>
                <DefaultButton onClick={() => saveSelectionB(selectedPoints)}>
                    update selection B ({selectionB.length})
                </DefaultButton>
            </div>
            <div className='flex flex-wrap items-center mb-2 justify-start'>
                <DefaultButton onClick={() => compareClusterExplanation(selectionA, selectionB)}>
                    compute
                </DefaultButton>
            </div>

            {explanation.features.length > 0 &&
                <p className="text-sm text-gray-500 w-fit mt-3 text-left" >
                    Features that are most different between&nbsp;<a
                        onClick={() => selectPointsScatterplot(explanation.selectionA, scatterplot)}
                        className="underline cursor-pointer">
                        selection A
                    </a>
                    &nbsp;and&nbsp;
                    <a
                        onClick={() => selectPointsScatterplot(explanation.selectionB, scatterplot)}
                        className="underline cursor-pointer">
                        selection B
                    </a>.
                </p>
            }
            {explanation.features.length > 0 && <p className="text-sm text-gray-500 w-fit mb-2" >
                Click on a feature to change the point colors.</p>}
            <div className="flex flex-col items-left justify-between">
                {explanation.features.map((feature) => {
                    return (
                        <span key={feature} onClick={() => pointColorOnChange(feature)}
                            className="text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                            {feature}
                        </span>
                    );
                })}
            </div>

        </span>
    )

}

export const Explanation = ({ exclus, selectedPoints, pointColorOnChange, pointColorOptions, scatterplot }) => {
    const [explanation, setExplanation] = useState({ "features": [], "method": null, "selectedPoints": [] });

    const explainabilityOptions = [
        { 'value': 'wasserstein', 'label': 'wasserstein' },
        { 'value': 'histogram', 'label': 'histogram' }
    ];
    const [explainabilityMethod, setExlainabilityMethod] = useState(explainabilityOptions[0].label);

    const getExplanation = (selectedPoints) => {
        explainCluster(selectedPoints, explainabilityMethod).then((explanation) => {
            const newExplanation = { "features": explanation, "method": explainabilityMethod, "selectedPoints": selectedPoints };
            setExplanation(newExplanation);
        });
    };


    useEffect(() => {
        setExplanation({ "features": [], "method": null, "selectedPoints": [] });
    }, [pointColorOptions]);

    return (
        <div>
            {exclus != null &&
                <span><h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-3" >InfoClus</h4>
                    <div className="flex flex-col items-left my-2 justify-between select-text">
                        {exclus != null && exclus.map(element => {
                            return <><p className="text-md  font-bold text-gray-900 text-left my-2">{element[0]}</p>
                                <p className="text-sm text-gray-900 text-left my-2 inline align-left">
                                    {
                                        element[1].map(subelement => {
                                            return <>{subelement}, </>
                                        })
                                    }</p></>
                        })}
                    </div>
                </span>
            }


            <h4 className="text-md font-large leading-6 text-gray-900 w-fit" >
                Cluster vs. Rest
            </h4>
            <p className="text-sm text-gray-500 text-left my-1">
                Select a group of points to compute features that stand out compared to the rest of the data.
            </p>

            {/* Explainability method */}
            <div className="flex flex-col items-left my-2 justify-between">
                <label className="text-sm text-gray-500 w-fit min-w-fit" >method</label>
                <ReactSelect
                    options={explainabilityOptions}
                    selected={explainabilityMethod}
                    onChange={(newValue) => { setExlainabilityMethod(newValue) }}
                    menuPlacement={'top'}
                />
            </div>

            <div className='flex flex-wrap items-center mb-2 justify-start'>
                <DefaultButton onClick={() => getExplanation(selectedPoints)}>
                    compute
                </DefaultButton>
            </div>
            {explanation.features.length > 0 && <p className="text-sm text-gray-500 w-fit mt-3" >
                Method {explanation.method},&nbsp;
                <a
                    onClick={() => selectPointsScatterplot(explanation.selectedPoints, scatterplot)}
                    className="underline cursor-pointer">
                    selected points&nbsp;
                </a>
                {explanation.selectedPoints.length}</p>
            }
            {explanation.features.length > 0 && <p className="text-sm text-gray-500 w-fit mb-2" >
                Click on a feature to change the point colors.</p>}
            <div className="flex flex-col items-left justify-between">
                {explanation.features.map((feature) => {
                    return (
                        <span key={feature} onClick={() => pointColorOnChange(feature)}
                            className="text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                            {feature}
                        </span>
                    );
                })}
            </div>

            <ComparingClusters
                selectedPoints={selectedPoints}
                pointColorOnChange={pointColorOnChange}
                pointColorOptions={pointColorOptions}
                scatterplot={scatterplot} />

        </div >
    );
}


