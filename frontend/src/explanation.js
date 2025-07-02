import { useEffect, useState } from "react";
import { explainCluster, compareClusters } from "./api";
import { DefaultButton } from "./buttons";
import { ReactSelect } from "./utils";
import { Histogram } from "./histogram";

const selectPointsScatterplot = (selectedPoints, scatterplot) => {
    scatterplot.deselect();
    scatterplot.select(selectedPoints);
}

const upArrow = <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor" className="size-3">
    <path strokeLinecap="round" strokeLinejoin="round" d="M4.5 10.5 12 3m0 0 7.5 7.5M12 3v18" />
</svg>
const downArrow = <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor" className="size-3">
    <path strokeLinecap="round" strokeLinejoin="round" d="M19.5 13.5 12 21m0 0-7.5-7.5M12 21V3" />
</svg>

const FeatureTags = ({ features, selectedPointColor, onFeatureClick }) => {
    return (
        <div className="flex flex-wrap gap-2 items-start justify-start">
            {features.map((feature) => {
                let name = feature[0]
                const isHigher = feature[1]
                const isSelected = name === selectedPointColor

                return (
                    <span
                        key={name}
                        onClick={() => onFeatureClick(name)}
                        className={`
                            inline-flex items-center gap-1 px-3 py-1 rounded-full text-sm cursor-pointer
                            transition-all duration-200 hover:shadow-md
                            bg-gray-200/20 text-gray-500 border border-gray-200
                            ${isSelected ? 'font-bold ring-2 ring-blue-400' : ''}
                        `}
                    >
                        {name}
                        {isHigher ? upArrow : downArrow}
                    </span>
                );
            })}
        </div>
    );
};

const ClusterVsRest = ({ selectedPoints, pointColorOnChange, pointColorOptions, pointColors, scatterplot, selectedPointColor, setIsLoading }) => {
    const [explanation, setExplanation] = useState({ "features": [], "method": null, "selectedPoints": [] });

    const explainabilityOptions = [
        { 'value': 'wasserstein', 'label': 'wasserstein' },
    ];
    const [explainabilityMethod, setExlainabilityMethod] = useState(explainabilityOptions[0].label);

    const getExplanation = (selectedPoints) => {
        setIsLoading(true);
        explainCluster(selectedPoints, explainabilityMethod).then((explanation) => {
            var features = explanation["features"].map((feature, i) => { return ([feature, explanation["higher_mean"][i]]) });

            const newExplanation = {
                "features": features,
                "method": explainabilityMethod,
                "selectedPoints": selectedPoints
            };
            setExplanation(newExplanation);

            setIsLoading(false);
        });
    };

    useEffect(() => {
        setExplanation({ "features": [], "method": null, "selectedPoints": [] });
    }, [pointColorOptions]);

    return (
        <div>
            <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-4" >
                Cluster vs. Rest
            </h4>
            <p className="text-sm text-gray-500 text-left my-1">
                Select a group of points to compute features that stand out compared to the rest of the data.
            </p>
            {/* Selection Controls */}
            <div className='flex flex-wrap items-center gap-2 mb-2'>
                <span className="min-w-[100px] h-fit select-none justify-center items-center 
                    rounded-md bg-gray-300 text-gray-500 mr-2 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
                    border-0 drop-shadow-sm border border-gray-600">
                    current selection ({selectedPoints.length})
                </span>
                <DefaultButton
                    onClick={() => getExplanation(selectedPoints)}
                    disabled={selectedPoints.length <= 2}
                >
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
            <FeatureTags features={explanation.features} selectedPointColor={selectedPointColor} onFeatureClick={pointColorOnChange} />

            {pointColors["values"].length > 0 &&
                explanation.features.some(feature => feature[0] === selectedPointColor) &&
                <Histogram
                    featureValues={pointColors["values"]}
                    xlabel={selectedPointColor}
                    selectedPoints={explanation.selectedPoints}
                    selectedGroupName="selected"
                    otherGroupName="not selected" />
            }
        </div>
    );
};


export const ComparingClusters = ({ selectedPoints, pointColorOnChange, pointColorOptions, pointColors, scatterplot, selectedPointColor, setIsLoading }) => {
    const [explanation, setExplanation] = useState({
        "features": [],
        "selectionA": [],
        "selectionB": []
    });

    const [selectionA, setSelectionA] = useState([]);
    const [selectionB, setSelectionB] = useState([]);
    const [compareMode, setCompareMode] = useState('rest'); // 'rest' or 'manual'

    const saveSelectionA = (selectedPoints) => {
        const filteredBPoints = selectionB.filter(point => !selectedPoints.includes(point));
        setSelectionA(selectedPoints);
        setSelectionB(filteredBPoints);
        setExplanation({
            "features": [],
            "selectionA": [],
            "selectionB": []
        });
    };

    const saveSelectionB = (selectedPoints) => {
        const filteredSelection = selectedPoints.filter(point => !selectionA.includes(point));
        setSelectionB(filteredSelection);
        setExplanation({
            "features": [],
            "selectionA": [],
            "selectionB": []
        });
    };

    useEffect(() => {
        setExplanation({
            "features": [],
            "selectionA": [],
            "selectionB": []
        });
        setSelectionA([]);
        setSelectionB([]);
        setCompareMode('rest');
    }, [pointColorOptions]);

    const compareClusterExplanation = (selectionA, selectionB) => {
        if (selectionA.length <= 2) {
            return;
        }

        let finalSelectionB = selectionB;
        // if comparing to rest or selectionB is empty, use all remaining points
        if (compareMode === 'rest' || selectionB.length == 0) {
            const restOfPoints = pointColors["values"].map((_, i) => i).filter(i => !selectionA.includes(i));
            setSelectionB(restOfPoints);
            finalSelectionB = restOfPoints;
        }

        setIsLoading(true);
        compareClusters(selectionA, finalSelectionB).then((explanation) => {
            var features = explanation["features"].map((feature, i) => { return ([feature, explanation["higher_mean"][i]]) });

            const newExplanation = {
                "features": features,
                "selectionA": selectionA,
                "selectionB": finalSelectionB
            };
            setExplanation(newExplanation);
            setIsLoading(false);
        });
    }

    return (
        <span>
            <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-3" >
                Cluster vs. Cluster
            </h4>           
            <p className="text-sm text-gray-500 text-left my-1">
                Select two groups of points to compute features that stand out in one group compared to the other.
                If points are in both groups, they will be removed from selection B.
            </p>
            
            {/* Compare Mode Toggle */}
            <div className="flex items-center gap-4 mb-3">
                <span className="text-sm text-gray-700 font-medium">Comparison mode</span>
                <label className="flex items-center gap-2 cursor-pointer">
                    <input
                        type="radio"
                        name="compareMode"
                        value="rest"
                        checked={compareMode === 'rest'}
                        onChange={(e) => setCompareMode(e.target.value)}
                        className="text-purple-600"
                    />
                    <span className="text-sm text-gray-600">rest of data</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer">
                    <input
                        type="radio"
                        name="compareMode"
                        value="manual"
                        checked={compareMode === 'manual'}
                        onChange={(e) => setCompareMode(e.target.value)}
                        className="text-purple-600"
                    />
                    <span className="text-sm text-gray-600">manual</span>
                </label>
            </div>
            
            <div className='flex justify-between items-center mb-2'>
                <span
                    className={`
                            flex flex-col items-start gap-1 px-3 py-1 text-sm rounded-sm
                            transition-all duration-200 
                            text-gray-500 border-l-4 border-[#800080] w-[47%]'
                        `}
                >
                    <div>selection A ({selectionA.length})</div>
                    <div className="flex justify-between w-full gap-1">
                        <button className='text-sm underline' onClick={() => saveSelectionA(selectedPoints)}>update</button>
                        <button className='text-sm underline' onClick={() => selectPointsScatterplot(selectionA, scatterplot)}>show</button>
                        <button className='text-sm underline' onClick={() => saveSelectionA([])}>clear</button>
                    </div>
                </span>
                {compareMode === 'manual' && (
                    <span
                        className={`
                                flex flex-col items-start gap-1 px-3 py-1 text-sm w-[47%] rounded-sm
                                transition-all duration-200 
                                text-gray-500 border-l-4 border-gray-400
                            `}
                    >
                        <div>selection B ({selectionB.length})</div>
                        <div className="flex justify-between w-full gap-1">
                            <button className='text-sm underline' onClick={() => saveSelectionB(selectedPoints)}>update</button>
                            <button className='text-sm underline' onClick={() => selectPointsScatterplot(selectionB, scatterplot)}>show</button>
                            <button className='text-sm underline' onClick={() => saveSelectionB([])}>clear</button>
                        </div>
                    </span>
                )}
            </div>

            {/* Selection Controls */}
            <div className='flex flex-wrap items-center gap-2 mb-2'>
                <DefaultButton
                    onClick={() => compareClusterExplanation(selectionA, selectionB)}
                    disabled={selectionA.length <= 2 || (compareMode === "manual" && selectionB.length <= 2)}
                >
                    compute
                </DefaultButton>
            </div>

            {/* Results Summary */}
            {explanation.features.length > 0 && (
                <p className="text-sm text-gray-500 w-fit mt-3 mb-2 text-left">
                    Arrows denote if the mean of a feature is higher in selection A. 
                    Click on a feature to change the point colors.
                </p>
            )}
            <FeatureTags features={explanation.features} selectedPointColor={selectedPointColor} onFeatureClick={pointColorOnChange} />

            {pointColors["values"] && pointColors["values"].length > 0 &&
                explanation.features.some(feature => feature[0] === selectedPointColor) &&
                (() => {
                    // Create filtered dataset with only selection A and B values
                    const combinedIndices = [...explanation.selectionA, ...explanation.selectionB];
                    const filteredValues = combinedIndices.map(index => pointColors["values"][index]);
                    const mappedSelectionA = Array.from({ length: explanation.selectionA.length }, (_, i) => i);

                    return (
                        <Histogram
                            featureValues={filteredValues}
                            xlabel={selectedPointColor}
                            selectedPoints={mappedSelectionA}
                            selectedGroupName="selection A"
                            selectedGroupColor="#800080"
                            otherGroupName="selection B" />
                    );
                })()
            }

        </span>
    )
}

export const Explanation = ({ exclus, selectedPoints, pointColorOnChange, pointColorOptions, pointColors, scatterplot, selectedPointColor }) => {
    const [isLoading, setIsLoading] = useState(false);

    return (
        <div>
            <p className="text-sm text-gray-500 text-left my-1">
                Rank features according to the Wasserstein distance between the histograms of
                the selected points and the rest of the data (cluster vs. rest) or between two groups of points (cluster vs. cluster).
            </p>

{/*             <ClusterVsRest
                selectedPoints={selectedPoints}
                pointColorOnChange={pointColorOnChange}
                pointColorOptions={pointColorOptions}
                pointColors={pointColors}
                scatterplot={scatterplot}
                selectedPointColor={selectedPointColor}
                setIsLoading={setIsLoading}
            /> */}

            <ComparingClusters
                selectedPoints={selectedPoints}
                pointColorOnChange={pointColorOnChange}
                pointColorOptions={pointColorOptions}
                pointColors={pointColors}
                scatterplot={scatterplot}
                selectedPointColor={selectedPointColor}
                setIsLoading={setIsLoading}
            />

            {isLoading && (
                <div
                    role="status" className="fixed top-0 left-0 w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center flex-col"
                >
                    <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600 m-4" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                    <span className="text-gray-500 text-sm">Computing features...</span>
                </div>
            )}
        </div >
    );
}


