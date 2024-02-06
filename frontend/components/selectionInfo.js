import { useState, useEffect } from "react";


// https://tailwindcomponents.com/component/profile-information-card-horizon-ui-tailwind
export default function SelectionInfo(props) {
    const { scatterplot, numPoints, datasetInfo } = props;

    const [numSelected, setNumSelected] = useState(0);
    useEffect(() => {
        if (scatterplot != null) {
            console.log("subscribing to select action");
            scatterplot.subscribe('select', onSelect);
            scatterplot.subscribe('deselect', onDeselect);
        } else {
            console.log("Info: scatterplot is null");
        }
    }, [scatterplot]);

    const onSelect = (points) => {
        if (points["points"].length > 0) {
            setNumSelected(points["points"].length)
        }
    }

    const onDeselect = () => {
        let cview = scatterplot.get('cameraView');
        setNumSelected(0)
        scatterplot.set({
            opacityBy: 'density',
            cameraView: cview,
        })
        scatterplot.refresh()
    }

    return (
        <div className="select-none absolute bottom-2 rounded-lg right-2 w-fit bg-white/80 p-1 flex text-sm text-gray-500"
        >
            {numSelected > 0 &&
                <p>{numSelected} selected</p>
            }
            <abbr title="hold shift to draw a selection, Shift+Ctrl adds to the current selection, double click on the background removes the selection">
                <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-5 h-5 mx-1">
                    <path fillRule="evenodd" d="M2.25 12c0-5.385 4.365-9.75 9.75-9.75s9.75 4.365 9.75 9.75-4.365 9.75-9.75 9.75S2.25 17.385 2.25 12zm8.706-1.442c1.146-.573 2.437.463 2.126 1.706l-.709 2.836.042-.02a.75.75 0 01.67 1.34l-.04.022c-1.147.573-2.438-.463-2.127-1.706l.71-2.836-.042.02a.75.75 0 11-.671-1.34l.041-.022zM12 9a.75.75 0 100-1.5.75.75 0 000 1.5z" clipRule="evenodd" />
                </svg>
            </abbr>
            <p>{numPoints} points, {datasetInfo}</p>
        </div>
    )

}