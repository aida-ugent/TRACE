import { useEffect, useState } from "react";
import { explainCluster, compareClusters } from "./api";
import { DefaultButton } from "./buttons";
import { ReactSelect } from "./utils";


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


export const ComparingClusters = ({ selectedPoints, pointColorOnChange, pointColorOptions, scatterplot, setIsLoading, selectedPointColor }) => {
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
        setIsLoading(true);
        compareClusters(selectionA, selectionB).then((explanation) => {
            var features = explanation["features"].map((feature, i) => { return ([feature, explanation["higher_mean"][i]]) });

            const newExplanation = {
                "features": features,
                "selectionA": selectionA,
                "selectionB": selectionB
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
                    </a>. Arrows denote if the mean of a feature is higher in selection A.
                </p>
            }
            {explanation.features.length > 0 && <p className="text-sm text-gray-500 w-fit mb-2" >
                Click on a feature to change the point colors.</p>}
            <div className="flex flex-col items-left justify-between">
                {explanation.features.map((feature) => {
                    let name = feature[0]
                    if (name == selectedPointColor) {
                        return (
                            <span key={name} onClick={() => pointColorOnChange(name)}
                                className="flex flex-row items-center text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                                <p className="font-bold">{name}</p>{feature[1] ? upArrow : downArrow}
                            </span>
                        );
                    } else {
                        return (
                            <span key={name} onClick={() => pointColorOnChange(name)}
                                className="flex flex-row items-center text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                                {name}{feature[1] ? upArrow : downArrow}
                            </span>
                        );
                    }
                })}
            </div>

        </span>
    )

}

export const Explanation = ({ exclus, selectedPoints, pointColorOnChange, pointColorOptions, scatterplot, selectedPointColor }) => {
    const [explanation, setExplanation] = useState({ "features": [], "method": null, "selectedPoints": [] });
    const [isLoading, setIsLoading] = useState(false);

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
            <p className="text-sm text-gray-500 text-left my-1">
                Rank features according to the Wasserstein distance between the histograms of
                the selected points and the rest of the data (cluster vs. rest) or between two groups of points (cluster vs. cluster).
            </p>

            <h4 className="text-md font-large leading-6 text-gray-900 w-fit mt-4" >
                Cluster vs. Rest
            </h4>
            <p className="text-sm text-gray-500 text-left my-1">
                Select a group of points to compute features that stand out compared to the rest of the data.
            </p>

            {/* Explainability method */}
            {/* <div className="flex flex-col items-left my-2 justify-between">
                <label className="text-sm text-gray-500 w-fit min-w-fit" >method</label>
                <ReactSelect
                    options={explainabilityOptions}
                    selected={explainabilityMethod}
                    onChange={(newValue) => { setExlainabilityMethod(newValue) }}
                    menuPlacement={'top'}
                />
            </div> */}

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
                    let name = feature[0]
                    if (name == selectedPointColor) {
                        return (
                            <span key={name} onClick={() => pointColorOnChange(name)}
                                className="flex flex-row items-center text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                                <p className="font-bold">{name}</p>{feature[1] ? upArrow : downArrow}
                            </span>
                        );
                    } else {
                        return (
                            <span key={name} onClick={() => pointColorOnChange(name)}
                                className="flex flex-row items-center text-sm text-gray-500 cursor-pointer w-fit min-w-fit">
                                {name}{feature[1] ? upArrow : downArrow}
                            </span>
                        );
                    }
                })}
            </div>

            <ComparingClusters
                selectedPoints={selectedPoints}
                pointColorOnChange={pointColorOnChange}
                pointColorOptions={pointColorOptions}
                scatterplot={scatterplot}
                setIsLoading={setIsLoading}
                selectedPointColor={selectedPointColor} />
            <div
                role="status" className="fixed top-0 left-0 w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center flex-col"
                style={{ visibility: isLoading ? "visible" : "hidden" }}
            >
                <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600 m-4" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                <span className="text-gray-500 text-sm">Computing features...</span>
            </div>
        </div >
    );
}


