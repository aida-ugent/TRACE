import { useState, useEffect, useCallback, useRef } from "react";
import { Tooltip } from 'react-tooltip';


// https://tailwindcomponents.com/component/profile-information-card-horizon-ui-tailwind
export default function SelectionInfo(props) {
    const { scatterplot, numPoints, datasetInfo, opacityByDensity, opacity } = props;

    const [numSelected, setNumSelected] = useState(0);
    
    // Use refs to store current values to avoid recreating callbacks
    const opacityRef = useRef(opacity);
    const opacityByDensityRef = useRef(opacityByDensity);
    
    // Update refs when props change
    opacityRef.current = opacity;
    opacityByDensityRef.current = opacityByDensity;
    
    const onSelect = useCallback((points) => {
        if (points["points"].length > 0) {
            setNumSelected(points["points"].length)
        }
    }, []);

    const onDeselect = useCallback(() => {
        setNumSelected(0)
        // only reset opacity if we showed the high dimensional neighbors
        if (scatterplot.get('opacityBy') === 'valueW') {            
            if (opacityByDensityRef.current) {
                scatterplot.set({
                    opacityBy: 'density',
                })
            } else {
                scatterplot.set({
                    opacityBy: null,
                    opacity: opacityRef.current['value'],
                })
            }
        }
    }, [scatterplot]); 

    useEffect(() => {
        if (scatterplot != null) {
            scatterplot.subscribe('select', onSelect);
            scatterplot.subscribe('deselect', onDeselect);
            
            // Cleanup function to unsubscribe when component unmounts
            return () => {
                scatterplot.unsubscribe('select', onSelect);
                scatterplot.unsubscribe('deselect', onDeselect);
            };
        } else {
            console.log("Info: scatterplot is null");
        }
    }, [scatterplot]);

    return (
        <div className="select-none absolute bottom-2 right-2 w-fit bg-white/80 p-1 flex text-sm"
        >
            {numSelected > 0 &&
                <p>{numSelected} selected</p>
            }
            <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor"
                className="size-5 cursor-pointer mr-1"
                data-tooltip-id="info-tooltip">
                <path strokeLinecap="round" strokeLinejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
            </svg>
            <Tooltip id="info-tooltip" className='max-w-[400px] text-sm text-left'>
                Hold shift to draw a selection, Shift+Ctrl adds to the current selection, double click on the background removes the selection.
            </Tooltip>
            <p>{numPoints} points, {datasetInfo}</p>
        </div>
    )

}